#!/usr/bin/env python3
"""Check the public source allowlist, empty study templates and release archives."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
IGNORED = {".git", "__pycache__", ".pytest_cache", ".ruff_cache", "build", "dist"}
FORBIDDEN_SUFFIXES = {".db", ".sqlite", ".duckdb", ".parquet", ".h5", ".sdf", ".mol2", ".pdb", ".pdbqt", ".cif", ".mmcif", ".xlsx", ".csv", ".tsv", ".pt", ".pth", ".ckpt", ".onnx", ".pkl", ".log", ".pem", ".key"}
PATTERNS = [
    re.compile(r"/(?:Users|home)/" + r"[A-Za-z][A-Za-z0-9_.-]*/"),
    re.compile(r"100\." + r"(?:6[4-9]|[7-9][0-9]|1[01][0-9]|12[0-7])\.[0-9]+\.[0-9]+(?!/10)"),
    re.compile(r"https://[a-z0-9-]+\." + r"[a-z0-9-]+\.chatgpt\.site"),
    re.compile(r"SkinPipe-(?:Collaborator|Research)-" + r"0\.2\.13rc1"),
    re.compile(r"(?:sk-proj-|sk-)" + r"[A-Za-z0-9_-]{35,}"),
    re.compile(r"-----BEGIN " + r"(?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(r"31" + r"개 (?:분석|연구|표적)|original" + r"31|routine200-" + r"qualification"),
]
SCREENSHOT_MANIFEST = "docs/images/manifest.json"


def screenshot_manifest(entries: dict[str, bytes]) -> tuple[dict[str, str], list[str]]:
    """Only individually reviewed, hash-pinned tutorial images are exempted."""
    matches = [name for name in entries if name == SCREENSHOT_MANIFEST or name.endswith("/" + SCREENSHOT_MANIFEST)]
    if not matches:
        return {}, []
    if len(matches) != 1:
        return {}, ["Package has more than one screenshot manifest"]
    manifest_name = matches[0]
    prefix = manifest_name[:-len(SCREENSHOT_MANIFEST)]
    approved, errors = {}, []
    try:
        value = json.loads(entries[manifest_name])
        assert value["scope"] == "public-tutorial-screenshots"
        for row in value["images"]:
            filename, digest = row["file"], row["sha256"]
            assert re.fullmatch(r"[0-9]{2}-[a-z0-9-]+\.(?:jpg|png|svg)", filename)
            assert re.fullmatch(r"[0-9a-f]{64}", digest)
            name = prefix + "docs/images/" + filename
            assert name not in approved
            if name not in entries:
                errors.append(f"Reviewed screenshot is missing: {name}")
            approved[name] = digest
    except (AssertionError, KeyError, TypeError, ValueError):
        return {}, ["Invalid reviewed screenshot manifest"]
    return approved, errors


def inspect_bytes(name: str, data: bytes, approved_images: dict[str, str] | None = None) -> list[str]:
    errors = []
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name:
        errors.append(f"Unsafe file path: {name}")
    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        errors.append(f"Research/runtime file type: {name}")
    if approved_images and name in approved_images:
        if hashlib.sha256(data).hexdigest() != approved_images[name]:
            errors.append(f"Reviewed screenshot hash changed: {name}")
        if path.suffix == ".jpg":
            if not (data.startswith(b"\xff\xd8\xff") and data.endswith(b"\xff\xd9")):
                errors.append(f"Reviewed screenshot is not a JPEG: {name}")
            return errors
        if path.suffix == ".png":
            if not (data.startswith(b"\x89PNG\r\n\x1a\n") and data.endswith(b"\x00\x00\x00\x00IEND\xaeB`\x82")):
                errors.append(f"Reviewed screenshot is not a PNG: {name}")
            return errors
        if path.suffix == ".svg":
            try:
                if ET.fromstring(data).tag != "{http://www.w3.org/2000/svg}svg":
                    raise ValueError("Not an SVG root")
            except (ET.ParseError, ValueError):
                errors.append(f"Reviewed graph is not an SVG: {name}")
            # Vector exports remain subject to all private-content checks below.
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return errors + [f"Unexpected binary file: {name}"]
    for pattern in PATTERNS:
        if pattern.search(text):
            errors.append(f"Private/credential content pattern: {name}")
    if path.name in {"skin_targets.yaml", "cosmetic_actives.yaml"} and "config" in path.parts:
        field = "markers" if path.name == "skin_targets.yaml" else "actives"
        if text.strip() != f"version: 1\n{field}: []":
            errors.append(f"Nonempty study configuration: {name}")
    return errors


def check_templates(root: Path) -> list[str]:
    errors = []
    config = root / "skinpipe/config"
    # These are empty platform templates, not a published target or result set.
    for filename, field in (("skin_targets.yaml", "markers"), ("cosmetic_actives.yaml", "actives")):
        if (config / filename).read_text().strip() != f"version: 1\n{field}: []":
            errors.append(f"Study configuration must remain empty: {filename}")
    for filename, fields in (("direction_evidence.json", ["evidence"]), ("curated_literature.json", ["sources", "claims"])):
        value = json.loads((config / filename).read_text())
        if any(value.get(field) != [] for field in fields):
            errors.append(f"Study evidence must remain empty: {filename}")
    value = json.loads((config / "structure_weight_overrides.json").read_text())
    if value.get("nodes") != {} or value.get("automatic_application") is not False:
        errors.append("Study-specific ranking overrides are not public platform defaults")
    return errors


def check_tree(root: Path) -> tuple[list[str], int]:
    manifest = root / "PUBLIC_FILES.txt"
    allowed = set(manifest.read_text().splitlines())
    files = []
    errors = []
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if any(part in IGNORED or part.endswith(".egg-info") for part in relative.parts):
            continue
        if path.is_symlink():
            errors.append(f"Symlink is not a reviewed source file: {relative}")
        elif path.is_file():
            files.append(path)
            if relative.as_posix() not in allowed:
                errors.append(f"File is outside the reviewed allowlist: {relative}")
    entries = {path.relative_to(root).as_posix(): path.read_bytes() for path in files}
    approved_images, image_errors = screenshot_manifest(entries)
    errors.extend(image_errors)
    for name, data in entries.items():
        errors.extend(inspect_bytes(name, data, approved_images))
    actual = {path.relative_to(root).as_posix() for path in files}
    errors.extend(f"Allowlisted file is missing: {name}" for name in sorted(allowed - actual))
    errors.extend(check_templates(root))
    return errors, len(files)


def check_package(path: Path) -> tuple[list[str], int]:
    errors = []
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            entries = [(name, archive.read(name)) for name in archive.namelist() if not name.endswith("/")]
    else:
        with tarfile.open(path) as archive:
            entries = []
            for member in archive.getmembers():
                if member.isdir():
                    continue
                if not member.isfile():
                    errors.append(f"Package contains link or special file: {member.name}")
                    continue
                stream = archive.extractfile(member)
                assert stream is not None
                entries.append((member.name, stream.read()))
    approved_images, image_errors = screenshot_manifest(dict(entries))
    errors.extend(image_errors)
    for name, data in entries:
        errors.extend(inspect_bytes(name, data, approved_images))
    return errors, len(entries)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--package", type=Path)
    args = parser.parse_args()
    errors, count = check_package(args.package) if args.package else check_tree(args.root.resolve())
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Public source check passed: {count} files; no research/runtime files or configured study inputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
