from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_wham_window_provenance.py"
spec = importlib.util.spec_from_file_location("wham_provenance", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["wham_provenance"] = module
spec.loader.exec_module(module)


def payload(times: list[float], coordinates: list[float]) -> bytes:
    return "".join(
        f"{time:.6f} {coordinate:.8f}\n"
        for time, coordinate in zip(times, coordinates)
    ).encode("utf-8")


def archive_path(system: str, route: str, center: float) -> str:
    token = str(center).replace(".", "-")
    return (
        f"WHAM_PMF/histogram-collection-{system}-{route}/"
        f"{token}-test-dist_vs_t.dat"
    )


class WhamWindowProvenanceTests(unittest.TestCase):
    def test_parse_series_and_summary(self) -> None:
        raw = b"# comment\n0.0 1.0\n1.0 2.0\n2.0 3.0\n"
        times, coordinates = module.parse_series(raw, "x.dat")
        self.assertEqual(times, [0.0, 1.0, 2.0])
        self.assertEqual(coordinates, [1.0, 2.0, 3.0])
        summary, retained = module.series_summary(
            raw, "x.dat", "glucose", "US-1", 1.2
        )
        self.assertEqual(summary["sample_count"], 3)
        self.assertEqual(summary["duration"], 2.0)
        self.assertTrue(summary["time_step"]["all_equal_within_tolerance"])
        self.assertIsNone(retained)

    def test_rejects_nonmonotonic_time(self) -> None:
        with self.assertRaises(module.ProvenanceError):
            module.parse_series(b"0 1\n0 2\n", "x.dat")

    def test_series_comparisons(self) -> None:
        left = ([0.0, 1.0, 2.0], [1.0, 2.0, 3.0])
        same = ([0.0, 1.0, 2.0], [1.0, 2.0, 3.0])
        shifted = ([0.0, 1.0, 2.0], [2.0, 3.0, 4.0])
        exact = module.compare_series(left, same)
        self.assertTrue(exact["exact_time_grid"])
        self.assertTrue(exact["exact_coordinate_series"])
        different = module.compare_series(left, shifted)
        self.assertFalse(different["exact_coordinate_series"])
        self.assertAlmostEqual(
            different["coordinate_mean_absolute_difference"], 1.0
        )
        self.assertAlmostEqual(different["coordinate_pearson"] or 0.0, 1.0)

    def test_prefix_comparison(self) -> None:
        short = ([0.0, 1.0], [3.0, 4.0])
        long = ([0.0, 1.0, 2.0], [9.0, 8.0, 7.0])
        result = module.prefix_comparison(short, long)
        self.assertTrue(result["time_grid_is_exact_prefix"])
        self.assertFalse(result["coordinate_series_is_exact_prefix"])
        self.assertAlmostEqual(result["sample_fraction"], 2 / 3)
        self.assertAlmostEqual(result["duration_fraction"], 0.5)

    def test_interpolation_comparison(self) -> None:
        lower = ([0.0, 1.0], [1.0, 3.0])
        upper = ([0.0, 1.0], [3.0, 5.0])
        middle = ([0.0, 1.0], [2.0, 4.0])
        result = module.interpolation_comparison(lower, middle, upper)
        self.assertTrue(result["same_time_grid"])
        self.assertTrue(result["samplewise_arithmetic_mean_exact"])
        self.assertEqual(
            result[
                "samplewise_arithmetic_mean_maximum_absolute_difference"
            ],
            0.0,
        )

    def test_full_synthetic_audit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "WHAM_PMF.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as handle:
                for route_number in range(1, 4):
                    route = f"US-{route_number}"
                    times = [float(index) for index in range(5)]
                    route_offset = route_number * 0.001
                    handle.writestr(
                        archive_path("galactose", route, 0.0),
                        payload(times, [0.1 * route_number + value for value in times]),
                    )
                    handle.writestr(
                        archive_path("galactose", route, 0.3),
                        payload(times, [0.2 * route_number + value for value in times]),
                    )
                    handle.writestr(
                        archive_path("galactose", route, 2.1),
                        payload(
                            times,
                            [
                                2.1 + route_offset + 0.01 * value
                                for value in times
                            ],
                        ),
                    )
                    handle.writestr(
                        archive_path("galactose", route, 2.25),
                        payload(
                            times,
                            [
                                2.26 + route_offset + 0.02 * value
                                for value in times
                            ],
                        ),
                    )
                    handle.writestr(
                        archive_path("galactose", route, 2.4),
                        payload(
                            times,
                            [
                                2.4 + route_offset + 0.03 * value
                                for value in times
                            ],
                        ),
                    )

                lower_times = [0.0, 1.0, 2.0, 3.0, 4.0]
                short_times = lower_times[:3]
                handle.writestr(
                    archive_path("fructofuranose", "US-2", 12.3),
                    payload(lower_times, [12.3, 12.4, 12.5, 12.6, 12.7]),
                )
                handle.writestr(
                    archive_path("fructofuranose", "US-2", 12.6),
                    payload(short_times, [12.55, 12.65, 12.75]),
                )
                handle.writestr(
                    archive_path("fructofuranose", "US-2", 12.9),
                    payload(lower_times, [12.9, 13.0, 13.1, 13.2, 13.3]),
                )

            original_size = module.ARCHIVE_SIZE_BYTES
            original_md5 = module.ARCHIVE_MD5
            original_sha256 = module.ARCHIVE_SHA256
            try:
                module.ARCHIVE_SIZE_BYTES = archive.stat().st_size
                module.ARCHIVE_MD5 = module.digest_file(archive, "md5")
                module.ARCHIVE_SHA256 = module.digest_file(archive, "sha256")
                output = root / "output.json"
                report = module.audit(archive, output)
            finally:
                module.ARCHIVE_SIZE_BYTES = original_size
                module.ARCHIVE_MD5 = original_md5
                module.ARCHIVE_SHA256 = original_sha256

            self.assertEqual(report["window_file_count"], 18)
            self.assertTrue(report["decision"]["all_window_files_parsed"])
            self.assertFalse(
                report["decision"]
                ["galactose_0_0_and_0_3_wham_series_are_exact_duplicates_in_any_route"]
            )
            self.assertTrue(
                report["decision"]
                ["galactose_2_25_windows_are_distinct_deposited_series_under_tested_exact_checks"]
            )
            self.assertFalse(
                report["decision"]
                ["galactose_2_25_windows_are_exact_duplicates_or_exact_samplewise_interpolations"]
            )
            self.assertFalse(
                report["decision"]["fructofuranose_US2_12_6_has_205_samples"]
            )
            self.assertTrue(
                report["decision"]
                ["fructofuranose_US2_12_6_time_grid_is_exact_prefix_of_both_neighbours"]
            )
            self.assertFalse(
                report["decision"]
                ["fructofuranose_US2_12_6_coordinate_series_is_exact_prefix_of_either_neighbour"]
            )
            self.assertEqual(json.loads(output.read_text()), report)

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.ProvenanceError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
