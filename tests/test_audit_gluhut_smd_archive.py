from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_smd_archive.py"
spec = importlib.util.spec_from_file_location("smd_inventory", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["smd_inventory"] = module
spec.loader.exec_module(module)


class SMDArchiveInventoryTests(unittest.TestCase):
    def test_classification(self):
        self.assertIn("topology", module.classify(".prmtop"))
        self.assertIn("coordinate", module.classify(".pdb"))
        self.assertIn("parameter", module.classify(".frcmod"))
        self.assertIn("protocol_or_documentation", module.classify(".mdin"))
        self.assertIn("trajectory", module.classify(".nc"))
        self.assertEqual(module.classify(".unknown"), ["other"])

    def test_unsafe_member_rejected(self):
        for name in ("../escape", "/absolute", "a\\b"):
            with self.assertRaises(module.InventoryError):
                module.safe_member(name)

    def test_synthetic_archive_inventory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive_path = root / "fixture.zip"
            with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("system/host.prmtop", "%VERSION\n")
                archive.writestr("system/host.rst7", "title\n")
                archive.writestr("system/host.frcmod", "MASS\n")
                archive.writestr("system/run.mdin", "&cntrl\n/")
                archive.writestr("system/traj.nc", b"binary trajectory")
            original_size = module.ARCHIVE_SIZE_BYTES
            original_md5 = module.ARCHIVE_MD5
            try:
                module.ARCHIVE_SIZE_BYTES = archive_path.stat().st_size
                module.ARCHIVE_MD5 = module.digest_file(archive_path, "md5")
                output = root / "inventory.json"
                report = module.audit(archive_path, output)
            finally:
                module.ARCHIVE_SIZE_BYTES = original_size
                module.ARCHIVE_MD5 = original_md5
            self.assertTrue(report["exact_reusable_input_bundle_present"])
            self.assertTrue(
                report["decision"]["can_begin_exact_parent_parameter_audit_from_this_archive_alone"]
            )
            self.assertEqual(report["classification_counts"]["trajectory"], 1)
            self.assertTrue(output.exists())
            loaded = json.loads(output.read_text())
            self.assertEqual(loaded["archive"]["file_count"], 5)

    def test_duplicate_member_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive_path = root / "duplicate.zip"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr("same.txt", "first")
                archive.writestr("same.txt", "second")
            original_size = module.ARCHIVE_SIZE_BYTES
            original_md5 = module.ARCHIVE_MD5
            try:
                module.ARCHIVE_SIZE_BYTES = archive_path.stat().st_size
                module.ARCHIVE_MD5 = module.digest_file(archive_path, "md5")
                with self.assertRaises(module.InventoryError):
                    module.audit(archive_path, root / "output.json")
            finally:
                module.ARCHIVE_SIZE_BYTES = original_size
                module.ARCHIVE_MD5 = original_md5


if __name__ == "__main__":
    unittest.main()
