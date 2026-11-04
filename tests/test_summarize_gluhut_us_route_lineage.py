from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "summarize_gluhut_us_route_lineage.py"
spec = importlib.util.spec_from_file_location("us_lineage", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["us_lineage"] = module
spec.loader.exec_module(module)


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def route_path(system: str, route: str, center_tenths: int) -> str:
    whole, fraction = divmod(center_tenths, 10)
    if system == "galactose":
        run = route[3]
        label = route[4]
        directory = f"run {run} ({label})"
    else:
        run = route[3]
        directory = f"run{run}"
    return (
        f"US_Trajectories/{system}/{directory}/"
        f"prod-{whole}-{fraction}-100-frames.binpos"
    )


def synthetic_records() -> tuple[dict, dict]:
    routes = [
        *(('glucose', f'run{run}') for run in range(1, 4)),
        *(('fructopyranose', f'run{run}') for run in range(1, 4)),
        *(('fructofuranose', f'run{run}') for run in range(1, 4)),
        ('galactose', 'run1b'),
        ('galactose', 'run2b'),
        ('galactose', 'run3b'),
        ('galactose', 'run4a'),
        ('galactose', 'run5a'),
    ]
    members = []
    for system, route in routes:
        for center_tenths in module.EXPECTED_CENTRES_TENTHS:
            signature_route = 'run4a' if route == 'run5a' else route
            content_label = f"{system}:{signature_route}:{center_tenths}"
            if (
                system == 'galactose'
                and route in {'run1b', 'run2b'}
                and center_tenths == 3
            ):
                content_label = f"{system}:{route}:0"
            members.append(
                {
                    'path': route_path(system, route, center_tenths),
                    'kind': 'file',
                    'extension': '.binpos',
                    'classifications': ['trajectory', 'scripps_binpos'],
                    'sha256': digest(content_label),
                    'uncompressed_size_bytes': (
                        module.BINPOS_HEADER_BYTES
                        + 100 * module.BINPOS_RECORD_BYTES
                    ),
                }
            )
    inventory = {
        'archive': {
            'filename': 'US_Trajectories.zip',
            'size_bytes': 123756824,
            'md5': 'c460e6a1ee97ac531fc284d8f07e381e',
            'sha256': 'a' * 64,
        },
        'members': members,
    }

    wham_route_specs = [
        *(('glucose', f'US-{run}') for run in range(1, 4)),
        *(('fructopyranose', f'US-{run}') for run in range(1, 4)),
        *(('fructofuranose', f'US-{run}') for run in range(1, 4)),
        *(('galactose', f'US-{run}') for run in range(1, 4)),
        ('galactose', 'US-routeA'),
    ]
    wham = {
        'archive': {
            'filename': 'WHAM_PMF.zip',
            'size_bytes': 41560830,
            'md5': 'dbc6723b8adcb85deb8230f63b6a51da',
            'sha256': 'b' * 64,
        },
        'routes': [
            {
                'system': system,
                'route': route,
                'window_count': 51,
                'centres_A': list(module.EXPECTED_CENTRES_A),
                'frames_minimum_window': 15001,
                'frames_median_window': 15001,
                'frames_maximum_window': 15001,
            }
            for system, route in wham_route_specs
        ],
    }
    return inventory, wham


class UmbrellaRouteLineageTests(unittest.TestCase):
    def test_path_and_center_parser(self) -> None:
        match = module.TRAJECTORY_RE.fullmatch(
            'US_Trajectories/galactose/run 4 (a)/'
            'prod-12-6-100-frames.binpos'
        )
        assert match
        self.assertEqual(match.group('system'), 'galactose')
        self.assertEqual(match.group('run'), '4')
        self.assertEqual(match.group('label'), 'a')
        self.assertEqual(
            module.parse_center_tenths(
                match.group('whole'), match.group('fraction')
            ),
            126,
        )

    def test_wham_mapping(self) -> None:
        self.assertEqual(
            module.wham_target('glucose', 2, None),
            ('glucose', 'US-2'),
        )
        self.assertEqual(
            module.wham_target('galactose', 3, 'b'),
            ('galactose', 'US-3'),
        )
        self.assertEqual(
            module.wham_target('galactose', 5, 'a'),
            ('galactose', 'US-routeA'),
        )
        with self.assertRaises(module.LineageError):
            module.wham_target('galactose', 4, 'b')

    def test_frame_count_from_size(self) -> None:
        size = module.BINPOS_HEADER_BYTES + 100 * module.BINPOS_RECORD_BYTES
        self.assertEqual(module.infer_binpos_frames(size), 100)
        with self.assertRaises(module.LineageError):
            module.infer_binpos_frames(size + 1)

    def test_full_synthetic_lineage(self) -> None:
        inventory, wham = synthetic_records()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inventory_path = root / 'inventory.json'
            wham_path = root / 'wham.json'
            output_path = root / 'lineage.json'
            inventory_path.write_text(json.dumps(inventory), encoding='utf-8')
            wham_path.write_text(json.dumps(wham), encoding='utf-8')

            report = module.run(inventory_path, wham_path, output_path)

            self.assertEqual(report['trajectory_file_count'], 714)
            self.assertEqual(report['raw_route_directory_count'], 14)
            self.assertEqual(
                report['unique_complete_route_content_signature_count'],
                13,
            )
            self.assertEqual(report['deposited_wham_route_count'], 13)
            decision = report['decision']
            self.assertTrue(decision['all_714_trajectory_files_parsed'])
            self.assertTrue(
                decision['expected_fourteen_route_directories_present']
            )
            self.assertTrue(
                decision['all_routes_have_51_windows_0_to_15_by_0_3A']
            )
            self.assertTrue(
                decision['all_trajectory_files_have_100_recorded_frames']
            )
            self.assertTrue(
                decision['all_route_centres_match_mapped_wham_centres']
            )
            self.assertTrue(
                decision[
                    'fourteen_route_directories_reduce_to_thirteen_unique_content_signatures'
                ]
            )
            self.assertTrue(
                decision[
                    'thirteen_unique_route_signatures_match_thirteen_wham_routes'
                ]
            )
            self.assertTrue(
                decision[
                    'galactose_run4a_and_run5a_are_exact_duplicate_routes'
                ]
            )
            self.assertTrue(
                decision[
                    'galactose_routeA_has_two_archive_copies_but_one_wham_route'
                ]
            )
            self.assertTrue(
                decision['within_route_duplicate_window_content_detected']
            )
            self.assertFalse(
                decision[
                    'umbrella_archive_contains_full_wham_histogram_frame_series'
                ]
            )
            self.assertTrue(
                decision[
                    'median_histogram_to_visualization_frame_ratio_is_approximately_150'
                ]
            )
            self.assertFalse(
                decision['histogram_generation_reproducible_from_us_archive']
            )
            self.assertFalse(
                decision['pmf_recalculation_from_us_trajectories_authorized']
            )

            duplicate_routes = report['duplicate_complete_route_groups']
            self.assertEqual(len(duplicate_routes), 1)
            self.assertEqual(
                duplicate_routes[0]['route_keys'],
                ['galactose/run4a', 'galactose/run5a'],
            )
            within = report['within_route_duplicate_groups_by_route']
            self.assertEqual(set(within), {'galactose/run1b', 'galactose/run2b'})
            for groups in within.values():
                self.assertEqual(len(groups), 1)
                self.assertEqual(groups[0]['centres_A'], [0.0, 0.3])
            self.assertEqual(json.loads(output_path.read_text()), report)


if __name__ == '__main__':
    unittest.main()
