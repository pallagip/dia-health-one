#!/usr/bin/env python3
"""Correct public author-contact status after a rate-limited send attempt.

The umbrella provenance addendum is stored as a Gmail draft in the existing
Bristol thread. It has not been delivered. This updater makes that distinction
explicit and appends an idempotent delivery-status section to the public record.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

STATUS_VALUE = "draft_created_after_gmail_rate_limit"
DELIVERY_SECTION = """

## Delivery status

The umbrella-specific provenance addendum was **not sent** on 4 November 2026 because Gmail returned an outbound rate-limit error. The complete message was saved as a draft in the existing Wayne Powell / Anthony Davis thread.

Earlier source-parameter, SMD atom-selection and endpoint-provenance requests were sent successfully. The umbrella addendum remains queued and must not be represented as delivered until Gmail confirms a successful send.
"""


class CorrectionError(RuntimeError):
    pass


def correct_status(status: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(status))
    external = result.setdefault("external_execution", {})
    external["author_umbrella_lineage_clarification_sent"] = False
    external["author_umbrella_lineage_clarification_draft_created"] = True
    external["author_umbrella_lineage_clarification_status"] = STATUS_VALUE
    return result


def correct_record(text: str) -> str:
    if "## Delivery status" in text:
        if "was **not sent**" not in text:
            raise CorrectionError(
                "existing delivery-status section does not contain the expected boundary"
            )
        return text
    return text.rstrip() + DELIVERY_SECTION + "\n"


def run(status_path: Path, record_path: Path) -> dict[str, Any]:
    status = json.loads(status_path.read_text(encoding="utf-8"))
    corrected = correct_status(status)
    status_path.write_text(
        json.dumps(corrected, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    record_path.write_text(
        correct_record(record_path.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    return corrected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--status-json",
        type=Path,
        default=Path("research/status/status.json"),
    )
    parser.add_argument(
        "--record",
        type=Path,
        default=Path("research/GLUHUT_UMBRELLA_LINEAGE_AUDIT.md"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    status = run(args.status_json, args.record)
    print(
        json.dumps(
            {
                "schema": status.get("schema"),
                "external_execution": status["external_execution"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
