"""Synthetic tests: no stored research inputs, candidates or predictions."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import subprocess
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
        value["images"][0]["file"] = "01-example.pdb"
        self.assertTrue(guard.screenshot_manifest({
            guard.SCREENSHOT_MANIFEST: json.dumps(value).encode(), "docs/images/01-example.pdb": b"data"})[1])

    def test_png_and_svg_are_hash_bound_and_text_still_checked(self):
        import hashlib
        png = b"\x89PNG\r\n\x1a\nreviewed-tutorial\x00\x00\x00\x00IEND\xaeB`\x82"
        svg = b'<svg xmlns="http://www.w3.org/2000/svg"><title>Tutorial</title></svg>'
        for filename, data in (("21-example.png", png), ("22-example.svg", svg)):
            name = "docs/images/" + filename
            digest = hashlib.sha256(data).hexdigest()
            value = {"scope": "public-tutorial-screenshots", "images": [
                {"file": filename, "sha256": digest}]}
            approved, errors = guard.screenshot_manifest({
                guard.SCREENSHOT_MANIFEST: json.dumps(value).encode(), name: data})
            self.assertEqual(errors, [])
            self.assertEqual(guard.inspect_bytes(name, data, approved), [])
            self.assertTrue(guard.inspect_bytes(name, data + b"changed", approved))
            self.assertTrue(guard.inspect_bytes(name, b"wrong format", {name: hashlib.sha256(b"wrong format").hexdigest()}))
        self.assertTrue(guard.inspect_bytes("docs/images/23-unreviewed.png", png))
        private_svg = svg.replace(b"Tutorial", ("/Users/" + "researcher/private").encode())
        name = "docs/images/22-example.svg"
        self.assertTrue(guard.inspect_bytes(name, private_svg, {name: hashlib.sha256(private_svg).hexdigest()}))

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

    def test_package_rejects_duplicate_zip_entries_before_name_lookup(self):
        import warnings
        import zipfile
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "duplicate.whl"
            with warnings.catch_warnings(), zipfile.ZipFile(archive, "w") as output:
                warnings.simplefilter("ignore", UserWarning)
                output.writestr("skinpipe/code.py", ("sk-proj-" + "a" * 40))
                output.writestr("skinpipe/code.py", "pass")
            errors, _ = guard.check_package(archive)
            self.assertTrue(any("Duplicate" in error for error in errors), errors)

    def test_package_rejects_zip_links_and_special_files(self):
        import stat
        import zipfile
        with tempfile.TemporaryDirectory() as temporary:
            for mode in (stat.S_IFLNK, stat.S_IFIFO, stat.S_IFCHR):
                with self.subTest(mode=mode):
                    archive = Path(temporary) / "special.whl"
                    info = zipfile.ZipInfo("skinpipe/code.py")
                    info.create_system = 3
                    info.external_attr = (mode | 0o600) << 16
                    with zipfile.ZipFile(archive, "w") as output:
                        output.writestr(info, "synthetic-target")
                    errors, _ = guard.check_package(archive)
                    self.assertTrue(any("special file" in error for error in errors), errors)

    def test_package_directory_names_are_checked(self):
        import zipfile
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "directory.whl"
            with zipfile.ZipFile(archive, "w") as output:
                output.writestr("../outside/", "")
                output.writestr("skinpipe/code.py", "pass")
            errors, _ = guard.check_package(archive)
            self.assertTrue(any("Unsafe" in error for error in errors), errors)

    def test_package_rejects_nonempty_study_json_templates(self):
        import zipfile
        values = {
            "direction_evidence.json": {"evidence": [{"synthetic": "study"}]},
            "curated_literature.json": {"sources": ["synthetic"], "claims": []},
            "structure_weight_overrides.json": {"nodes": {"synthetic": 2}, "automatic_application": False},
        }
        with tempfile.TemporaryDirectory() as temporary:
            for name, value in values.items():
                with self.subTest(name=name):
                    archive = Path(temporary) / "study.whl"
                    with zipfile.ZipFile(archive, "w") as output:
                        output.writestr("skinpipe/_assets/config/" + name, json.dumps(value))
                    errors, _ = guard.check_package(archive)
                    self.assertTrue(any("Study" in error for error in errors), errors)

    def test_package_limits_are_checked_before_payload_read(self):
        import zipfile
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "large.whl"
            with zipfile.ZipFile(archive, "w") as output:
                output.writestr("skinpipe/code.py", "pass")
            with patch.object(guard, "MAX_PACKAGE_BYTES", 3, create=True), patch.object(zipfile.ZipFile, "read", side_effect=AssertionError("payload must not be read")):
                errors, _ = guard.check_package(archive)
            self.assertTrue(any("limit" in error for error in errors), errors)

    def test_tar_package_duplicates_and_file_parent_conflicts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            for names in (("source/code.py", "source/code.py"), ("source/code.py", "source/code.py/child")):
                with self.subTest(names=names):
                    archive = Path(temporary) / "ambiguous.tar.gz"
                    with tarfile.open(archive, "w:gz") as output:
                        for name in names:
                            info = tarfile.TarInfo(name)
                            info.size = 4
                            output.addfile(info, io.BytesIO(b"pass"))
                    errors, _ = guard.check_package(archive)
                    self.assertTrue(any("Duplicate" in error or "conflict" in error for error in errors), errors)

    def test_empty_templates_and_regular_tar_directories_remain_valid(self):
        values = {"direction_evidence.json": {"evidence": []}, "curated_literature.json": {"sources": [], "claims": []},
                  "structure_weight_overrides.json": {"nodes": {}, "automatic_application": False}}
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "source.tar.gz"
            with tarfile.open(archive, "w:gz") as output:
                directory = tarfile.TarInfo("source/config")
                directory.type = tarfile.DIRTYPE
                output.addfile(directory)
                for name, value in values.items():
                    data = json.dumps(value).encode()
                    info = tarfile.TarInfo("source/config/" + name)
                    info.size = len(data)
                    output.addfile(info, io.BytesIO(data))
            self.assertEqual(guard.check_package(archive), ([], 3))

    def test_corrupt_package_returns_a_failure_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "broken.tar.gz"
            archive.write_bytes(b"not-an-archive")
            errors, count = guard.check_package(archive)
            self.assertEqual(count, 0)
            self.assertTrue(any("Unreadable" in error for error in errors), errors)

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

    def test_radar_reference_tracks_skin_scoring_policy(self):
        import ast
        import re
        module = (ROOT / "web/public/metric-radar.js").read_text()
        reference = json.loads(re.search(r"export const SKIN_REFERENCE = (\{.*?\});", module, re.S).group(1))
        policy = ast.parse((ROOT / "skinpipe/skinpipe/filters/topical.py").read_text())
        constants = {}
        for node in policy.body:
            if isinstance(node, ast.Assign) and len(node.targets) == 1:
                target = node.targets[0]
                if isinstance(target, ast.Tuple) and all(isinstance(n, ast.Name) for n in target.elts):
                    constants.update(zip([n.id for n in target.elts], ast.literal_eval(node.value)))
        self.assertEqual(reference["molecular_weight"]["max"], constants["_MW_IDEAL"])
        self.assertEqual([reference["logP"]["min"], reference["logP"]["max"]], [constants["_LOGP_LO"], constants["_LOGP_HI"]])
        scoring = next(node for node in policy.body if isinstance(node, ast.FunctionDef) and node.name == "topical_score")
        for variable, field in [("s_hbd", "hydrogen_bond_donors"), ("s_tpsa", "tpsa")]:
            value = next(node.value for node in scoring.body if isinstance(node, ast.Assign) and node.targets[0].id == variable)
            self.assertEqual(ast.literal_eval(value.body), 10.0)
            self.assertEqual(reference[field]["max"], ast.literal_eval(value.test.comparators[0]))
        self.assertNotIn("QED", reference)
        self.assertNotIn("hydrogen_bond_acceptors", reference)


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

    def test_archive_duplicate_paths_are_rejected_before_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "duplicate.tar.gz"
            with tarfile.open(archive, "w:gz") as output:
                for name in ("source/code.py", "source/./code.py"):
                    info = tarfile.TarInfo(name)
                    info.size = 4
                    output.addfile(info, io.BytesIO(b"pass"))
            with self.assertRaisesRegex(ValueError, "Duplicate|canonical"):
                bootstrap.safe_unpack(archive, root / "output")
            self.assertFalse((root / "output").exists())

    def test_archive_file_parent_conflict_and_size_limit_preserve_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "source.tar.gz"
            with tarfile.open(archive, "w:gz") as output:
                for name in ("source/code.py", "source/code.py/child"):
                    info = tarfile.TarInfo(name)
                    info.size = 4
                    output.addfile(info, io.BytesIO(b"pass"))
            with self.assertRaisesRegex(ValueError, "conflict"):
                bootstrap.safe_unpack(archive, root / "output")
            with patch.object(bootstrap, "MAX_ARCHIVE_BYTES", 3), self.assertRaisesRegex(ValueError, "limit"):
                bootstrap.safe_unpack(archive, root / "output")
            self.assertFalse((root / "output").exists())

    def test_rapid_source_retries_preserve_each_interrupted_download(self):
        def fail_download(url, target):
            target.write_bytes(b"next interrupted download")
            raise OSError("synthetic interrupted download")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "source.tar.gz").write_bytes(b"first interrupted download")
            with patch.object(bootstrap.dt, "datetime") as clock, patch.object(bootstrap, "download", side_effect=fail_download), patch.object(bootstrap.urllib.request, "urlopen") as request:
                clock.now.return_value.strftime.return_value = "20000101000000"
                request.side_effect = lambda *args, **kwargs: io.BytesIO(json.dumps({"sha": "a" * 40}).encode())
                for _ in range(2):
                    with self.assertRaisesRegex(OSError, "synthetic"):
                        bootstrap.prepare_source(root, None, {})
            preserved = list(root.glob("interrupted-source-*.tar.gz"))
            self.assertEqual(len(preserved), 2)
            self.assertEqual({p.read_bytes() for p in preserved}, {b"first interrupted download", b"next interrupted download"})
            self.assertFalse((root / "source").exists())

    def prepared(self, temporary, ready):
        root = Path(temporary).resolve() / "app"
        root.mkdir()
        (root / "source/skinpipe").mkdir(parents=True)
        (root / "source/skinpipe/pyproject.toml").write_text("synthetic")
        (root / "artifacts").mkdir()
        (root / "artifacts/app.whl").write_bytes(b"synthetic wheel")
        state = {"schema": bootstrap.SCHEMA, "library": "skin-starter", "device": "auto", "source_ready": True, "wheel": "app.whl", "wheel_sha256": bootstrap.hashlib.sha256(b"synthetic wheel").hexdigest(), "catalog_ready": ready}
        bootstrap.write_json(root / "installation.json", state)
        (root / "runtime/app").mkdir(parents=True)
        for name in (
            "Prepare-Catalog.sh", "Catalog-Status.sh", "Start-SkinPipe.sh",
            "Stop-SkinPipe.sh", "app/Start-SkinPipe.sh", "app/Stop-SkinPipe.sh",
        ):
            launcher = root / "runtime" / name
            launcher.write_text("#!/usr/bin/env bash\nexit 0\n")
            launcher.chmod(0o700)
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

    @patch.object(bootstrap, "check_platform")
    @patch.object(bootstrap, "run")
    def test_incomplete_launchers_do_not_resume_or_claim_ready(self, run, platform):
        for ready in (False, True):
            for missing in (True, False):
                with self.subTest(ready=ready, missing=missing), tempfile.TemporaryDirectory() as temporary:
                    root, args = self.prepared(temporary, ready)
                    launcher = root / "runtime/app/Start-SkinPipe.sh"
                    if missing:
                        launcher.unlink()
                    else:
                        launcher.chmod(0o600)
                    with self.assertRaisesRegex(ValueError, "launchers are incomplete"):
                        bootstrap.install(args)
                    run.assert_not_called()
                    self.assertEqual(
                        json.loads((root / "installation.json").read_text())["catalog_ready"],
                        ready,
                    )

    def test_failed_launcher_permissions_do_not_publish_resume_entry(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run(
                ["bash", "-c", '''
source "$1"
task_root="$2"
task_library=skin-starter
task_device=cpu
chmod() { return 17; }
if task_write_prepare_catalog; then exit 99; fi
[[ ! -e "$task_root/Prepare-Catalog.sh" ]]
''', "launcher-test", str(ROOT / "scripts/install-from-sources-ubuntu.sh"), temporary],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((Path(temporary) / "Prepare-Catalog.sh.tmp").is_file())


if __name__ == "__main__":
    unittest.main()
