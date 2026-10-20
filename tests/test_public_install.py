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

    def test_screenshot_exception_is_exact_and_hash_bound(self):
        import hashlib
        name = "docs/images/01-example.jpg"
        data = b"\xff\xd8\xffreviewed-tutorial\xff\xd9"
        approved = {name: hashlib.sha256(data).hexdigest()}
        self.assertEqual(guard.inspect_bytes(name, data, approved), [])
        self.assertTrue(guard.inspect_bytes(name, data + b"changed", approved))
        self.assertTrue(guard.inspect_bytes("docs/images/02-unreviewed.jpg", data, approved))
        self.assertTrue(guard.inspect_bytes(name, data))

    def test_screenshot_manifest_cannot_authorize_other_files(self):
        value = {"scope": "public-tutorial-screenshots", "images": [
            {"file": "../research.jpg", "sha256": "0" * 64}]}
        self.assertTrue(guard.screenshot_manifest({guard.SCREENSHOT_MANIFEST: json.dumps(value).encode()})[1])
        value["images"][0]["file"] = "01-example.jpg"
        self.assertTrue(guard.screenshot_manifest({guard.SCREENSHOT_MANIFEST: json.dumps(value).encode()})[1])

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

    def test_metric_help_is_available_without_sorting(self):
        import re
        module = (ROOT / "web/public/column-help.js").read_text()
        payload = re.search(r"export const metricColumnHelp = (\{.*?\});", module, re.S)
        self.assertIsNotNone(payload)
        help_text = json.loads(payload.group(1))
        self.assertIn("↑", help_text["priority_score_v2"][2])
        self.assertIn("↓", help_text["original_analysis_rank"][2])
        self.assertIn("↓", help_text["review_dock_affinity"][2])
        self.assertIn("같은 수용체", help_text["review_dock_affinity"][2])
        self.assertIn("적정 범위", help_text["logp"][2])
        self.assertNotIn("↑", help_text["permeation"][2])
        self.assertIn('aria-hidden="true">?</span>', module)
        self.assertIn("'pointerover'", module)
        self.assertIn("'focusin'", module)
        self.assertIn("button.closest('dialog[open]')", module)
        self.assertIn("event.preventDefault()", module)
        # Help controls have no data-sort attribute and cannot invoke sorting.
        button = module[module.index("export function helpButton"):module.index("export function installColumnHelp")]
        self.assertNotIn("data-sort", button)
        entry = (ROOT / "web/public/local-research-app.js").read_text()
        self.assertIn("installColumnHelp();startLocalResearch()", entry)

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
