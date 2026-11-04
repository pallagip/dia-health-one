from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_us_archive.py"
spec = importlib.util.spec_from_file_location("us_inventory", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["us_inventory"] = module
spec.loader.exec_module(module)


class UmbrellaArchiveInventoryTests(unittest.TestCase):
    def test_system_route_and_window_inference(self) -> None:
        path = "US_Trajectories/glucose/US-2/window-12.6/traj.binpos"
        self.assertEqual(module.infer_system(path), "glucose")
        self.assertEqual(module.infer_route(path), "2")
        self.assertIn(12.6, module.infer_window_candidates(path))

    def test_binary_format_detection(self) -> None:
        self.assertEqual(
            module.detect_binary_format(b"fxyz" + b"\0" * 16, ".binpos"),
            "scripps_binpos",
        )
        self.assertEqual(
            module.detect_binary_format(b"CDF\x01" + b"\0" * 16, ".nc"),
            "netcdf_classic",
        )
        self.assertEqual(
            module.detect_binary_format(b"\x89HDF\r\n\x1a\n" + b"\0" * 16, ".nc"),
            "hdf5_or_netcdf4",
        )

    def test_rst_content_classification(self) -> None:
        restraint = b"&rst\n iat=-1,-1,\n igr1=1,2,\n igr2=3,4,\n/\n"
        self.assertIn("restraint_input", module.classify(".rst", restraint))
        self.assertIn(
            "protocol_or_documentation",
            module.classify(".rst", restraint),
        )
        self.assertIn(
            "ambiguous_rst_or_restart",
            module.classify(".rst", b"title\n2\n0 0 0 1 1 1\n"),
        )

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../escape", "/absolute", "a\\b"):
            with self.assertRaises(module.InventoryError):
                module.safe_member(value)

    def test_fixture_inventory_detects_duplicate_content(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive_path = root / "fixture.zip"
            with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr(
                    "US_Trajectories/glucose/US-1/window-1.5/a.binpos",
                    b"fxyz" + b"\0" * 12,
                )
                archive.writestr(
                    "US_Trajectories/glucose/US-1/window-1.5/b.binpos",
                    b"fxyz" + b"\0" * 12,
                )
                archive.writestr(
                    "US_Trajectories/glucose/US-1/run.mdin",
                    "&cntrl\n/\n",
                )
                archive.writestr(
                    "US_Trajectories/glucose/US-1/system.parm7",
                    "%VERSION\n",
                )

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

            self.assertTrue(report["archive"]["crc_verified"])
            self.assertEqual(report["archive"]["file_count"], 4)
            self.assertTrue(
                report["decision"]["exact_duplicate_file_content_detected"]
            )
            self.assertEqual(
                report["duplicate_file_content_groups"][0]["file_count"],
                2,
            )
            self.assertEqual(
                report["system_file_counts"],
                {"glucose": 4},
            )
            self.assertEqual(
                json.loads(output.read_text(encoding="utf-8")),
                report,
            )

    def test_duplicate_member_path_rejected(self) -> None:
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
