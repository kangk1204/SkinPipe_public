"""Resolve audited offline assets without depending on the caller's cwd."""
from pathlib import Path

REQUIRED_ASSETS = (
    "bridge/export_images.py",
    "public/run-offline.html",
    "public/run-offline.js",
    "public/run-offline-worker.js",
    "public/vendor/molstar/molstar.js",
    "public/vendor/molstar/molstar.css",
    "public/vendor/molstar/LICENSE",
    "public/vendor/molstar/NOTICE.json",
)


def export_asset_root():
    module = Path(__file__).resolve()
    packaged = module.parent / "_export_assets"
    root = packaged if packaged.is_dir() else module.parents[3] / "web"
    if any(not (root / name).is_file() for name in REQUIRED_ASSETS):
        raise ValueError("Offline export assets missing; reinstall from a complete SkinPipe wheel or checkout")
    return root


def implementation_files():
    module = Path(__file__).resolve()
    root = export_asset_root()
    files = {"skinpipe/skinpipe/platform/" + name: module.with_name(name)
             for name in ("exports.py", "pose_view.py", "export_assets.py")}
    files.update({"web/" + name: root / name for name in REQUIRED_ASSETS})
    return files
