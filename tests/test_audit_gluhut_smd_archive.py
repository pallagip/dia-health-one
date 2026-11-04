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
    def test_classification_uses_rst_content(self):
        self.assertIn("topology", module.classify(".prmtop"))
        self.assertIn("parameterized_topology", module.classify(".parm7"))
        self.assertIn("ambiguous_rst", module.classify(".rst"))
        self.assertIn(
            "restraint_input",
            module.classify(".rst", b"&rst\n iat=-1,-1,\n r1=1.0,\n/\n"),
        )
        self.assertIn(
            "protocol_or_documentation",
            module.classify(".rst", b"&rst\n iat=-1,-1,\n/\n"),
        )
        self.assertIn(
            "coordinate",
            module.classify(".rst", b"title\n2\n0 0 0 1 1 1\n"),
        )
        self.assertIn("coordinate", module.classify(".rst7"))
        self.assertIn("source_parameter", module.classify(".frcmod"))
        self.assertIn("protocol_or_documentation", module.classify(".mdin"))
        self.assertIn("trajectory", module.classify(".binpos"))
        self.assertIn("trajectory", module.classify(".nc"))
        self.assertEqual(module.classify(".unknown"), ["other"])

    def test_unsafe_member_rejected(self):
        for name in ("../escape", "/absolute", "a\\b"):
            with self.assertRaises(module.InventoryError):
                module.safe_member(name)

    def _audit_fixture(self, members: dict[str, str | bytes]) -> dict:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        archive_path = root / "fixture.zip"
        with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
            for path, payload in members.items():
                archive.writestr(path, payload)
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
        self.assertTrue(output.exists())
        self.assertEqual(json.loads(output.read_text()), report)
        return report

    def test_compiled_trajectory_bundle_is_distinct_from_source_provenance(self):
        report = self._audit_fixture(
            {
                "system/host.parm7": "%VERSION\n%FLAG CHARGE\n",
                "system/distance.rst": "&rst\n iat=-1,-1,\n r1=1.0,\n/\n",
                "system/run.in": "&cntrl\n/\n",
                "system/traj.binpos": b"binary trajectory",
            }
        )
        components = report["exact_input_components"]
        self.assertTrue(components["topology_present"])
        self.assertTrue(components["parameterized_topology_present"])
        self.assertFalse(components["coordinate_present"])
        self.assertTrue(components["restraint_input_present"])
        self.assertTrue(components["protocol_present"])
        self.assertTrue(components["trajectory_present"])
        self.assertFalse(components["source_parameter_files_present"])

        self.assertTrue(report["compiled_topology_trajectory_bundle_present"])
        self.assertFalse(report["standalone_coordinate_bundle_present"])
        self.assertFalse(report["parameter_generation_provenance_complete"])
        self.assertFalse(report["exact_reusable_input_bundle_present"])
        self.assertTrue(report["decision"]["can_begin_compiled_topology_audit"])
        self.assertTrue(
            report["decision"]["can_begin_compiled_topology_and_trajectory_audit"]
        )
        self.assertFalse(
            report["decision"]["can_begin_standalone_restart_coordinate_audit"]
        )
        self.assertFalse(
            report["decision"]
            ["can_reconstruct_parameter_generation_from_this_archive"]
        )
        self.assertEqual(report["classification_counts"]["trajectory"], 1)
        self.assertEqual(report["classification_counts"]["restraint_input"], 1)
        self.assertNotIn("coordinate", report["classification_counts"])

    def test_coordinate_and_source_parameters_complete_broader_bundle(self):
        report = self._audit_fixture(
            {
                "system/host.parm7": "%VERSION\n%FLAG CHARGE\n",
                "system/host.rst7": "title\n2\n0.0 0.0 0.0\n",
                "system/host.frcmod": "MASS\n",
                "system/host.mdin": "&cntrl\n/\n",
                "system/traj.nc": b"binary trajectory",
            }
        )
        self.assertTrue(report["compiled_topology_trajectory_bundle_present"])
        self.assertTrue(report["standalone_coordinate_bundle_present"])
        self.assertTrue(report["parameter_generation_provenance_complete"])
        self.assertTrue(report["exact_reusable_input_bundle_present"])
        self.assertTrue(
            report["decision"]["can_reconstruct_parameter_generation_from_this_archive"]
        )

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
