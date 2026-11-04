from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "correct_umbrella_contact_status.py"
spec = importlib.util.spec_from_file_location("contact_status", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["contact_status"] = module
spec.loader.exec_module(module)


class ContactStatusCorrectionTests(unittest.TestCase):
    def test_correct_status(self) -> None:
        status = {
            "schema": "dia-health-one.public-confidence-status.v10",
            "external_execution": {
                "author_umbrella_lineage_clarification_sent": True,
                "other": 1,
            },
        }
        result = module.correct_status(status)
        external = result["external_execution"]
        self.assertFalse(
            external["author_umbrella_lineage_clarification_sent"]
        )
        self.assertTrue(
            external["author_umbrella_lineage_clarification_draft_created"]
        )
        self.assertEqual(
            external["author_umbrella_lineage_clarification_status"],
            module.STATUS_VALUE,
        )
        self.assertEqual(external["other"], 1)
        self.assertTrue(status["external_execution"]["author_umbrella_lineage_clarification_sent"])

    def test_record_update_is_idempotent(self) -> None:
        original = "# Record\n\nContent.\n"
        first = module.correct_record(original)
        second = module.correct_record(first)
        self.assertEqual(first, second)
        self.assertIn("## Delivery status", first)
        self.assertIn("was **not sent**", first)
        self.assertIn("saved as a draft", first)

    def test_existing_incompatible_section_fails(self) -> None:
        with self.assertRaises(module.CorrectionError):
            module.correct_record(
                "# Record\n\n## Delivery status\n\nSent successfully.\n"
            )

    def test_full_file_update(self) -> None:
        status = {
            "schema": "dia-health-one.public-confidence-status.v10",
            "external_execution": {
                "author_umbrella_lineage_clarification_sent": True
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            status_path = root / "status.json"
            record_path = root / "record.md"
            status_path.write_text(json.dumps(status), encoding="utf-8")
            record_path.write_text("# Record\n", encoding="utf-8")
            result = module.run(status_path, record_path)
            self.assertFalse(
                result["external_execution"]
                ["author_umbrella_lineage_clarification_sent"]
            )
            self.assertEqual(json.loads(status_path.read_text()), result)
            self.assertIn("was **not sent**", record_path.read_text())



    def test_confirmed_delivery_preserved(self) -> None:
        status = {
            "external_execution": {
                "author_umbrella_lineage_clarification_sent": False,
                "author_umbrella_lineage_clarification_sent_message_id": module.SENT_MESSAGE_ID,
                "author_umbrella_lineage_clarification_sent_thread_id": module.SENT_THREAD_ID,
            },
        }
        first = module.correct_status(status)
        second = module.correct_status(first)
        self.assertEqual(first, second)
        external = first["external_execution"]
        self.assertTrue(external["author_umbrella_lineage_clarification_sent"])
        self.assertFalse(external["author_umbrella_lineage_clarification_draft_pending"])
        self.assertEqual(external["author_umbrella_lineage_clarification_status"], module.DELIVERED_VALUE)

    def test_sent_record_is_idempotent(self) -> None:
        first = module.correct_record("# Record\n", delivered=True)
        self.assertEqual(module.correct_record(first, delivered=True), first)
        self.assertIn("1a1a8f2999166f97", first)
        self.assertIn("Gmail send", first)
        with self.assertRaises(module.CorrectionError):
            module.correct_record(first, delivered=False)

    def test_unchanged_historical_note_is_replaced_after_gmail_confirmation(self) -> None:
        historical = module.correct_record("# Record\n", delivered=False)
        delivered = module.correct_record(historical, delivered=True)
        self.assertNotIn("was **not sent**", delivered)
        self.assertIn("**Delivered:**", delivered)

if __name__ == "__main__":
    unittest.main()
