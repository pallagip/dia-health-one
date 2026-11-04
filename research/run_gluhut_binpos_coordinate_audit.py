#!/usr/bin/env python3
"""Run the GluHUT coordinate audit with the local audited BINPOS reader."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import audit_gluhut_binpos_coordinates as core
from binpos_native_reader import READER_VERSION, read_binpos


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--candidates-json",
        type=Path,
        default=Path(
            "research/data/gluhut-fructofuranose-restraint-candidates.json"
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-binpos-coordinate-audit.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    core.read_binpos = read_binpos
    report = core.audit(args.archive, args.candidates_json, args.output_json)

    # Replace the historical field name emitted by the core audit with an exact
    # reader provenance record. No result values are changed.
    report["software"].pop("mdtraj", None)
    report["software"]["binpos_reader"] = READER_VERSION
    report["software"]["binpos_format_reference"] = (
        "Scripps/VMD binposplugin.c as historically vendored by MDTraj 1.9.7"
    )
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "software": report["software"],
                "trajectory_count": report["trajectory_count"],
                "total_frame_count": report["total_frame_count"],
                "system_summaries": report["system_summaries"],
                "fructofuranose_alternatives": report[
                    "fructofuranose_alternative_cross_run_summary"
                ],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
