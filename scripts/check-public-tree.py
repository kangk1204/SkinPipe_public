#!/usr/bin/env python3
"""Check the public source allowlist, empty study templates and release archives."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import stat
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
MAX_PACKAGE_BYTES = 300 * 1024**2
MAX_PACKAGE_MEMBERS = 20000


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


def path_errors(name: str) -> list[str]:
    path = PurePosixPath(name)
    if (not path.parts or path.is_absolute() or ".." in path.parts or "\\" in name
            or name.rstrip("/") != path.as_posix() or any(ord(c) < 32 or ord(c) == 127 for c in name)):
        return [f"Unsafe file path: {name!r}"]
    return []


def study_template_errors(name: str, text: str) -> list[str]:
    path = PurePosixPath(name)
    if "config" not in path.parts:
        return []
    if path.name in {"skin_targets.yaml", "cosmetic_actives.yaml"}:
        field = "markers" if path.name == "skin_targets.yaml" else "actives"
        if text.strip() != f"version: 1\n{field}: []":
            return [f"Nonempty study configuration: {name}"]
    fields = {"direction_evidence.json": ("evidence",), "curated_literature.json": ("sources", "claims")}
    if path.name in fields or path.name == "structure_weight_overrides.json":
        try:
            value = json.loads(text)
            valid = isinstance(value, dict)
            if path.name in fields:
                valid = valid and all(value.get(field) == [] for field in fields[path.name])
            else:
                valid = valid and value.get("nodes") == {} and value.get("automatic_application") is False
        except (ValueError, TypeError):
            valid = False
        if not valid:
            return [f"Study evidence/configuration must remain empty: {name}"]
    return []


def inspect_bytes(name: str, data: bytes, approved_images: dict[str, str] | None = None) -> list[str]:
    errors = path_errors(name)
    path = PurePosixPath(name)
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
    errors.extend(study_template_errors(name, text))
    return errors


def check_templates(root: Path) -> list[str]:
    required = ("skin_targets.yaml", "cosmetic_actives.yaml", "direction_evidence.json",
                "curated_literature.json", "structure_weight_overrides.json")
    errors = []
    for filename in required:
        path = root / "skinpipe/config" / filename
        if not path.is_file():
            errors.append(f"Required empty study template is missing: {filename}")
        else:
            errors.extend(study_template_errors(path.relative_to(root).as_posix(), path.read_text()))
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
    return list(dict.fromkeys(errors)), len(files)


def check_package(path: Path) -> tuple[list[str], int]:
    errors, entries, seen, files = [], [], set(), set()
    total = 0

    def admit(name: str, size: int, directory: bool, regular: bool) -> None:
        nonlocal total
        errors.extend(path_errors(name))
        identity = PurePosixPath(name).as_posix()
        if identity in seen:
            errors.append(f"Duplicate package member: {name!r}")
        seen.add(identity)
        if regular:
            files.add(PurePosixPath(name))
        if not directory and not regular:
            errors.append(f"Package contains link or special file: {name!r}")
        if size < 0 or (directory and size):
            errors.append(f"Invalid package member size: {name!r}")
        total += size

    try:
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                members = archive.infolist()
                if len(members) > MAX_PACKAGE_MEMBERS:
                    return ["Package exceeds the member limit"], 0
                for member in members:
                    mode = stat.S_IFMT(member.external_attr >> 16)
                    directory = member.is_dir()
                    valid_type = mode in ({0, stat.S_IFDIR} if directory else {0, stat.S_IFREG})
                    admit(member.filename, member.file_size, directory and valid_type, not directory and valid_type)
                    if member.flag_bits & 1:
                        errors.append(f"Encrypted package member: {member.filename!r}")
                if total > MAX_PACKAGE_BYTES:
                    errors.append("Package exceeds the uncompressed size limit")
                if any(parent in files for name in seen for parent in PurePosixPath(name).parents):
                    errors.append("Package file/directory paths conflict")
                if errors:
                    return errors, 0
                # Read the exact central-directory entry, never a lossy name lookup.
                entries = [(member.filename, archive.read(member)) for member in members if not member.is_dir()]
        else:
            with tarfile.open(path) as archive:
                members = []
                for member in archive:
                    members.append(member)
                    admit(member.name, member.size, member.isdir(), member.isfile())
                    if len(members) > MAX_PACKAGE_MEMBERS or total > MAX_PACKAGE_BYTES:
                        return errors + ["Package exceeds the size/member limit"], 0
                if any(parent in files for name in seen for parent in PurePosixPath(name).parents):
                    errors.append("Package file/directory paths conflict")
                if errors:
                    return errors, 0
                for member in members:
                    if member.isdir():
                        continue
                    stream = archive.extractfile(member)
                    if stream is None:
                        return [f"Unreadable package member: {member.name!r}"], 0
                    with stream:
                        entries.append((member.name, stream.read()))
    except (OSError, ValueError, RuntimeError, tarfile.TarError, zipfile.BadZipFile) as error:
        return [f"Unreadable package: {type(error).__name__}"], 0
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
