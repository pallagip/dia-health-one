#!/usr/bin/env python3
"""Reconcile Bristol umbrella-lineage author contact against recorded Gmail evidence.

A historic rate-limited draft must not overwrite a verified Gmail SENT receipt.
Without the exact confirmed message and thread IDs, keep the status pending.
This updater never sends email or interprets SENT as a read acknowledgement.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

STATUS_VALUE = "draft_created_after_gmail_rate_limit"
DELIVERED_VALUE = "sent_confirmed_by_gmail"
SENT_MESSAGE_ID = "1a1a8f2999166f97"
SENT_THREAD_ID = "1a1a7ce37fa8d5f5"
DELIVERY_SECTION_PENDING = """\n\n## Delivery status

The umbrella-specific provenance addendum was **not sent** on 4 November 2026 because Gmail returned an outbound rate-limit error. The complete message was saved as a draft in the existing Wayne Powell / Anthony Davis thread.

Earlier source-parameter, SMD atom-selection and endpoint-provenance requests were sent successfully. The umbrella addendum remains queued and must not be represented as delivered until Gmail confirms a successful send.
"""
DELIVERY_SECTION_SENT = """\n\n## Delivery status

**Delivered:** Gmail confirmed the umbrella-specific provenance addendum was sent on 4 November 2026 to Wayne Powell (CC: Anthony Davis) in the existing GluHUT reproducibility thread. The earlier 4 November Gmail rate limit was transient and the saved draft was successfully sent without editing it.

- Gmail sent-message ID: `1a1a8f2999166f97`.
- Gmail thread ID: `1a1a7ce37fa8d5f5`.
- Gmail-reported message timestamp: `2026-11-04T23:43:53` (provider response has no timezone suffix).
- Delivery evidence: Gmail returned message ID, matching thread ID and `SENT` label, and the sent message was independently read back.
- Original draft ID `r-943222025130158032` is no longer pending; the `draft_created` status field records its historical creation, not an unsent message.

This is evidence of a successful **Gmail send**, not of recipient delivery confirmation, reading, author acknowledgement, or a scientific answer. All source-lineage and parameterization gaps remain open until independently resolved.

"""

class CorrectionError(RuntimeError):
    pass


def correct_status(status: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(status))
    external = result.setdefault("external_execution", {})
    sent = (
        external.get("author_umbrella_lineage_clarification_sent_message_id")
        == SENT_MESSAGE_ID
        and external.get("author_umbrella_lineage_clarification_sent_thread_id")
        == SENT_THREAD_ID
    )
    external["author_umbrella_lineage_clarification_sent"] = sent
    external["author_umbrella_lineage_clarification_draft_created"] = True
    external["author_umbrella_lineage_clarification_draft_pending"] = not sent
    external["author_umbrella_lineage_clarification_status"] = (
        DELIVERED_VALUE if sent else STATUS_VALUE
    )
    return result


def correct_record(text: str, delivered: bool = False) -> str:
    target = DELIVERY_SECTION_SENT if delivered else DELIVERY_SECTION_PENDING
    expression = r"(?ms)^## Delivery status\n.*?(?=^## |\Z)"
    existing = re.search(expression, text)
    if existing:
        section = existing.group()
        if not delivered and "was **not sent**" not in section:
            raise CorrectionError("cannot regress an existing delivered section")
        prefix = text[:existing.start()].rstrip()
        suffix = text[existing.end():].lstrip("\n")
        return prefix + target + (suffix if suffix else "")
    return text.rstrip() + target


def run(status_path: Path, record_path: Path) -> dict[str, Any]:
    status = json.loads(status_path.read_text(encoding="utf-8"))
    corrected = correct_status(status)
    status_path.write_text(
        json.dumps(corrected, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    record_path.write_text(
        correct_record(
            record_path.read_text(encoding="utf-8"),
            delivered=corrected["external_execution"]["author_umbrella_lineage_clarification_sent"],
        ),
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
