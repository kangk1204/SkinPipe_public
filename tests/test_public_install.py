"""Synthetic tests: no stored research inputs, candidates or predictions."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bootstrap = load("skinpipe_install", ROOT / "install.py")
guard = load("skinpipe_public_guard", ROOT / "scripts/check-public-tree.py")


class PublicSourceTests(unittest.TestCase):
    def test_reviewed_tree(self):
        errors, count = guard.check_tree(ROOT)
        self.assertEqual(errors, [])
        self.assertGreater(count, 200)

    def test_guard_rejects_output_and_credentials(self):
        self.assertTrue(guard.inspect_bytes("result" + ".xlsx", b"data"))
        self.assertTrue(guard.inspect_bytes("file.txt", ("sk-proj-" + "a" * 40).encode()))
        self.assertTrue(guard.inspect_bytes("file.txt", ("/Users/" + "researcher/research").encode()))

    def test_package_guard(self):
        import zipfile
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "software.whl"
            with zipfile.ZipFile(archive, "w") as output:
                output.writestr("skinpipe/code.py", "print('ok')\n")
            self.assertEqual(guard.check_package(archive)[0], [])
            with zipfile.ZipFile(archive, "a") as output:
                output.writestr("research/data" + ".csv", "private")
            self.assertTrue(guard.check_package(archive)[0])

    def test_browser_imports_resolve(self):
        import re
        public = ROOT / "web/public"
        for module in public.glob("*.js"):
            for relative in re.findall(r"(?:from\s*|import\s*)['\"](\./[^'\"]+)['\"]", module.read_text()):
                self.assertTrue((module.parent / relative).is_file(), (module.name, relative))

    def test_source_versions_agree(self):
        init = (ROOT / "skinpipe/skinpipe/__init__.py").read_text()
        project = (ROOT / "skinpipe/pyproject.toml").read_text()
        self.assertIn('"0.3.0rc1"', init)
        self.assertIn('version = "0.3.0rc1"', project)


class InstallerTests(unittest.TestCase):
    def test_unknown_prefix_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "keep.txt").write_text("keep")
            with self.assertRaises(ValueError):
                bootstrap.checked_prefix(str(root))
            self.assertEqual((root / "keep.txt").read_text(), "keep")

    def test_symlink_and_git_paths_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "link").symlink_to(root, target_is_directory=True)
            with self.assertRaises(ValueError):
                bootstrap.checked_prefix(str(root / "link/app"))
            (root / ".git").mkdir()
            with self.assertRaises(ValueError):
                bootstrap.checked_prefix(str(root / "app"))

    def test_archive_traversal_and_links_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for index, name in enumerate(("source/../../outside", "source/link")):
                archive = root / f"bad-{index}.tar.gz"
                with tarfile.open(archive, "w:gz") as output:
                    info = tarfile.TarInfo(name)
                    if index:
                        info.type = tarfile.SYMTYPE
                        info.linkname = "/tmp"
                        output.addfile(info)
                    else:
                        info.size = 1
                        output.addfile(info, io.BytesIO(b"x"))
                with self.assertRaises(ValueError):
                    bootstrap.safe_unpack(archive, root / f"output-{index}")
                self.assertFalse((root / f"output-{index}").exists())

    def test_archive_regular_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "source.tar.gz"
            with tarfile.open(archive, "w:gz") as output:
                info = tarfile.TarInfo("source/skinpipe/code.py")
                info.size = 4
                output.addfile(info, io.BytesIO(b"pass"))
            bootstrap.safe_unpack(archive, root / "source")
            self.assertEqual((root / "source/skinpipe/code.py").read_text(), "pass")

    def prepared(self, temporary, ready):
        root = Path(temporary).resolve() / "app"
        root.mkdir()
        (root / "source/skinpipe").mkdir(parents=True)
        (root / "source/skinpipe/pyproject.toml").write_text("synthetic")
        (root / "artifacts").mkdir()
        (root / "artifacts/app.whl").write_bytes(b"synthetic wheel")
        state = {"schema": bootstrap.SCHEMA, "library": "skin-starter", "device": "auto", "source_ready": True, "wheel": "app.whl", "wheel_sha256": bootstrap.hashlib.sha256(b"synthetic wheel").hexdigest(), "catalog_ready": ready}
        bootstrap.write_json(root / "installation.json", state)
        (root / "runtime").mkdir()
        (root / "runtime/Prepare-Catalog.sh").write_text("synthetic")
        args = argparse.Namespace(prefix=str(root), library="skin-starter", device="auto", accept_source_terms=False, check=False, source=None)
        return root, args

    @patch.object(bootstrap, "check_platform")
    @patch.object(bootstrap, "run")
    def test_prepared_start_does_not_rebuild(self, run, platform):
        with tempfile.TemporaryDirectory() as temporary:
            root, args = self.prepared(temporary, True)
            self.assertEqual(bootstrap.install(args), root / "runtime")
            run.assert_not_called()

    @patch.object(bootstrap, "check_platform")
    @patch.object(bootstrap, "run")
    def test_interrupted_catalog_uses_same_launcher(self, run, platform):
        with tempfile.TemporaryDirectory() as temporary:
            root, args = self.prepared(temporary, False)
            bootstrap.install(args)
            run.assert_called_once_with(["bash", str(root / "runtime/Prepare-Catalog.sh")], root / "catalog-resume.log")
            self.assertTrue(json.loads((root / "installation.json").read_text())["catalog_ready"])

    @patch.object(bootstrap, "check_platform")
    def test_concurrent_install_rejected(self, platform):
        with tempfile.TemporaryDirectory() as temporary:
            root, args = self.prepared(temporary, False)
            with (root / ".installer.lock").open("a") as handle:
                bootstrap.fcntl.flock(handle, bootstrap.fcntl.LOCK_EX | bootstrap.fcntl.LOCK_NB)
                with self.assertRaisesRegex(ValueError, "Another installer"):
                    bootstrap.install(args)

    @patch.object(bootstrap, "check_platform")
    @patch.object(bootstrap, "run", side_effect=RuntimeError("synthetic failure"))
    def test_failed_catalog_not_marked_ready(self, run, platform):
        with tempfile.TemporaryDirectory() as temporary:
            root, args = self.prepared(temporary, False)
            with self.assertRaises(RuntimeError):
                bootstrap.install(args)
            self.assertFalse(json.loads((root / "installation.json").read_text())["catalog_ready"])


if __name__ == "__main__":
    unittest.main()
