#!/usr/bin/env python3
"""Install the public SkinPipe platform, prepare a catalog, and start its local UI."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request

try:
    import fcntl
except ImportError:  # The platform check provides the Windows instructions.
    fcntl = None

REPOSITORY = "kangk1204/SkinPipe_public"
SCHEMA = "skinpipe-public-install-v1"
MAX_ARCHIVE_BYTES = 300 * 1024**2


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.chmod(0o600)
    temporary.replace(path)


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def check_platform() -> None:
    values = {}
    if Path("/etc/os-release").is_file():
        for line in Path("/etc/os-release").read_text().splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                values[key] = value.strip('"')
    if (sys.platform != "linux" or platform.machine() != "x86_64"
            or values.get("ID") != "ubuntu" or values.get("VERSION_ID") != "24.04"):
        raise ValueError("Install on Ubuntu 24.04 x86_64. On macOS/Windows, open the Ubuntu server through an SSH tunnel; see the README.")
    if os.geteuid() == 0:
        raise ValueError("Run as your regular Ubuntu user. The installer requests sudo only for operating-system packages.")


def checked_prefix(value: str) -> Path:
    path = Path(os.path.abspath(os.path.expanduser(value)))
    if path == Path(path.anchor) or path == Path.home():
        raise ValueError("Use a dedicated installation folder, such as ~/SkinPipe.")
    for parent in (path, *path.parents):
        if parent.is_symlink():
            raise ValueError(f"Installation paths must not contain symbolic links: {parent}")
        if (parent / ".git").exists():
            raise ValueError("Keep the installation and research data outside a Git checkout.")
    if path.exists():
        marker = path / "installation.json"
        if not marker.is_file() or read_json(marker).get("schema") != SCHEMA:
            raise ValueError(f"Folder already exists and is not a recognized SkinPipe installation: {path}. Choose a new --prefix; nothing was overwritten.")
        if path.stat().st_uid != os.geteuid():
            raise ValueError("The installation must be owned by your current Ubuntu user.")
    return path


def safe_unpack(archive: Path, destination: Path) -> None:
    """Extract only regular GitHub source files; reject traversal and all links."""
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        total = sum(member.size for member in members)
        if total > MAX_ARCHIVE_BYTES or len(members) > 20000:
            raise ValueError("Source archive exceeds the installation limit.")
        roots = set()
        for member in members:
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts or "\\" in member.name or not name.parts:
                raise ValueError("Unsafe source archive path.")
            if not member.isfile() and not member.isdir():
                raise ValueError("Links and special files are not accepted in source archives.")
            roots.add(name.parts[0])
        if len(roots) != 1:
            raise ValueError("Expected one source archive root.")
        destination.mkdir()
        for member in members:
            relative = PurePosixPath(member.name).parts[1:]
            if not relative:
                continue
            target = destination.joinpath(*relative)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = handle.extractfile(member)
                if stream is None:
                    raise ValueError("Unreadable source archive entry.")
                with stream, target.open("xb") as output:
                    shutil.copyfileobj(stream, output)
                target.chmod(0o600)


def download(url: str, target: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "SkinPipe-public-installer"})
    with urllib.request.urlopen(request, timeout=60) as response, target.open("xb") as output:
        length = 0
        while chunk := response.read(1024 * 1024):
            length += len(chunk)
            if length > MAX_ARCHIVE_BYTES:
                raise ValueError("Source download exceeds the installation limit.")
            output.write(chunk)


def prepare_source(root: Path, local_source: Path | None, state: dict) -> Path:
    source = root / "source"
    if state.get("source_ready"):
        if not (source / "skinpipe/pyproject.toml").is_file():
            raise ValueError("Saved source is incomplete. Preserve this folder and use a new --prefix.")
        return source
    if source.exists():
        raise ValueError("An incomplete source folder is preserved. Use a new --prefix.")
    if local_source is not None:
        original = local_source.resolve(strict=True)
        subprocess.run([sys.executable, str(original / "scripts/check-public-tree.py"), "--root", str(original)], check=True)
        shutil.copytree(original, source, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc", ".pytest_cache", ".ruff_cache", "*.egg-info", "build", "dist"))
        state["source_commit"] = "local-reviewed-source"
    else:
        request = urllib.request.Request(f"https://api.github.com/repos/{REPOSITORY}/commits/main", headers={"User-Agent": "SkinPipe-public-installer"})
        with urllib.request.urlopen(request, timeout=30) as response:
            metadata = json.load(response)
            commit = metadata.get("sha") if isinstance(metadata, dict) else None
        if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
            raise ValueError("GitHub did not return a valid source revision.")
        archive = root / "source.tar.gz"
        if archive.exists():
            archive.replace(root / f"interrupted-source-{dt.datetime.now().strftime('%Y%m%d%H%M%S')}.tar.gz")
        download(f"https://codeload.github.com/{REPOSITORY}/tar.gz/{commit}", archive)
        safe_unpack(archive, source)
        state["source_commit"] = commit
        state["source_archive_sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
    state["source_ready"] = True
    write_json(root / "installation.json", state)
    return source


def run(command: list[str], log: Path) -> None:
    print(f"Log: {log}", flush=True)
    with log.open("a") as output:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=False)
        assert process.stdout is not None
        for line in process.stdout:
            output.write(line)
            output.flush()
            print(line, end="", flush=True)
        code = process.wait()
    if code:
        raise RuntimeError(f"Step exited with status {code}. See {log}. Saved files were preserved.")


def install(args: argparse.Namespace) -> Path:
    check_platform()
    root = checked_prefix(args.prefix)
    required = 300 if args.library == "full" else 20
    ancestor = next(parent for parent in (root, *root.parents) if parent.exists())
    if not root.exists() and shutil.disk_usage(ancestor).free < required * 1024**3:
        raise ValueError(f"At least {required} GiB free space is required for {args.library}.")
    if args.library == "full" and not args.accept_source_terms:
        raise ValueError("Read docs/catalog.en.md and the providers' terms first. Full-library preparation requires --accept-source-terms.")
    if args.check:
        print(f"Ubuntu and installation path checks passed. Library: {args.library}; required free space: {required} GiB; path: {root}")
        return root
    os.umask(0o077)
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    root.chmod(0o700)
    lock = (root / ".installer.lock").open("a+")
    try:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("Another installer is already using this folder. Follow its terminal; no second build was started.") from error
        marker = root / "installation.json"
        state = read_json(marker) if marker.exists() else {"schema": SCHEMA, "repository": REPOSITORY, "library": args.library, "device": args.device, "created_at": dt.datetime.now(dt.timezone.utc).isoformat()}
        if state.get("library") != args.library or state.get("device") != args.device:
            raise ValueError("Keep the same --library and --device when resuming. Use a new --prefix for a different installation.")
        write_json(marker, state)
        source = prepare_source(root, args.source, state)
        print("[1/4] Preparing Ubuntu packages and the SkinPipe application.", flush=True)
        artifacts = root / "artifacts"
        artifacts.mkdir(exist_ok=True)
        if not state.get("wheel"):
            run(["sudo", "apt-get", "update"], root / "system-install.log")
            run(["sudo", "apt-get", "install", "-y", "python3.12-venv"], root / "system-install.log")
            build = root / "build-env"
            if not build.exists():
                subprocess.run(["python3.12", "-m", "venv", str(build)], check=True)
            run([str(build / "bin/python"), "-m", "pip", "wheel", "--no-deps", "--wheel-dir", str(artifacts), str(source / "skinpipe")], root / "wheel-build.log")
            wheels = list(artifacts.glob("skinpipe-*.whl"))
            if len(wheels) != 1:
                raise ValueError("Expected exactly one SkinPipe wheel.")
            state["wheel"] = wheels[0].name
            state["wheel_sha256"] = hashlib.sha256(wheels[0].read_bytes()).hexdigest()
            write_json(marker, state)
        wheel = artifacts / state["wheel"]
        if hashlib.sha256(wheel.read_bytes()).hexdigest() != state["wheel_sha256"]:
            raise ValueError("The saved installation wheel has changed. Nothing was installed.")
        runtime = root / "runtime"
        if not state.get("catalog_ready"):
            print("[2/4] Downloading official data/model files. [3/4] Computing ADMET and building the search index.\nThe first preparation can take a long time. Full-library CPU calculations may take days. Terminal logs and Catalog-Status.sh show progress; the same command resumes saved catalog work.", flush=True)
            prepare = runtime / "Prepare-Catalog.sh"
            if prepare.is_file():
                run(["bash", str(prepare)], root / "catalog-resume.log")
            elif runtime.exists():
                raise ValueError(f"Dependency installation stopped before catalog preparation. Its files and logs remain in {runtime}. Use a new --prefix after fixing the logged error; no catalog was restarted.")
            else:
                run(["bash", str(source / "scripts/install-from-sources-ubuntu.sh"), str(wheel), str(runtime), args.library, args.device], root / "prepare.log")
            state["catalog_ready"] = True
            write_json(marker, state)
        print(f"[4/4] Ready. Browser: http://127.0.0.1:4317\nInstallation: {root}\nStart: bash '{runtime}/Start-SkinPipe.sh'\nStop: bash '{runtime}/Stop-SkinPipe.sh'", flush=True)
        return runtime
    finally:
        lock.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="~/SkinPipe")
    parser.add_argument("--library", choices=("skin-starter", "full"), default="skin-starter")
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--accept-source-terms", action="store_true", help="Acknowledge the source terms reviewed for full-library preparation.")
    parser.add_argument("--source", type=Path, help="Use a reviewed local public checkout instead of downloading GitHub source.")
    parser.add_argument("--check", action="store_true", help="Check platform/path/capacity without installing or starting anything.")
    parser.add_argument("--no-start", action="store_true", help="Prepare and verify the installation without starting the browser server.")
    args = parser.parse_args(argv)
    try:
        runtime = install(args)
        if not args.check and not args.no_start:
            os.execvp("bash", ["bash", str(runtime / "Start-SkinPipe.sh")])
        return 0
    except KeyboardInterrupt:
        print("\nInstallation interrupted. Saved files remain in the installation folder. Re-run the same command to resume catalog preparation.", file=sys.stderr)
        return 130
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f"SkinPipe: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
