"""Ship the same audited offline assets in wheels and source distributions.

The web checkout is the source of truth. Build copies go only into the build
or release tree; an sdist contains its own copy and needs no sibling checkout.
"""
from pathlib import Path
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist

ROOT = Path(__file__).resolve().parent
ASSETS = (
    "bridge/export_images.py",
    "public/run-offline.html",
    "public/run-offline.js",
    "public/run-offline-worker.js",
    "public/vendor/molstar/molstar.js",
    "public/vendor/molstar/molstar.css",
    "public/vendor/molstar/LICENSE",
    "public/vendor/molstar/NOTICE.json",
)


def copy_export_assets(destination):
    packaged = ROOT / "skinpipe/platform/_export_assets"
    source = packaged if packaged.is_dir() else ROOT.parent / "web"
    missing = [name for name in ASSETS if not (source / name).is_file()]
    if missing:
        raise RuntimeError("Build requires the complete SkinPipe checkout or sdist; missing export assets: " + ", ".join(missing))
    names = set(ASSETS)
    # Personal installation ships the browser app, without campaign result data.
    if (source / "public/local-research.html").is_file():
        names.update(str(path.relative_to(source)) for path in (source / "public").rglob("*")
                     if path.is_file() and not path.is_symlink()
                     and path.suffix in {".html", ".js", ".css", ".svg", ".png", ".ico", ".json"}
                     and "data" not in path.relative_to(source / "public").parts)
    for name in sorted(names):
        target = Path(destination) / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)


class BuildPy(build_py):
    def run(self):
        # Restored sources may have older mtimes than ignored build/lib copies.
        # Copy every declared module from source rather than trusting timestamps.
        self.force = True
        super().run()
        copy_export_assets(Path(self.build_lib) / "skinpipe/platform/_export_assets")


class SourceDist(sdist):
    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        copy_export_assets(Path(base_dir) / "skinpipe/platform/_export_assets")


setup(cmdclass={"build_py": BuildPy, "sdist": SourceDist})
