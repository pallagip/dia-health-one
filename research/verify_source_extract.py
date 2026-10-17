#!/usr/bin/env python3
"""Check CSV values and source-cell references against the publisher XLSX."""

import csv
import hashlib
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXPECTED_SHA256 = "09a9b6aa4103f9d61af7df03886d59f16d6b09f69f2fe768260548061c7dc81d"
NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_DOC_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PACKAGE_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 research/verify_source_extract.py /path/to/source-data.xlsx")
    workbook_path = Path(sys.argv[1])
    digest = hashlib.sha256(workbook_path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise SystemExit(f"Source workbook SHA-256 mismatch: {digest}")

    with zipfile.ZipFile(workbook_path) as workbook:
        root = ET.fromstring(workbook.read("xl/workbook.xml"))
        sheet = root.find(f"{{{NS_MAIN}}}sheets/{{{NS_MAIN}}}sheet[@name='Figure 3d']")
        if sheet is None:
            raise SystemExit("Worksheet 'Figure 3d' was not found")
        rel_id = sheet.attrib[f"{{{NS_DOC_REL}}}id"]
        rels = ET.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        target = next(item.attrib["Target"] for item in rels.findall(f"{{{NS_PACKAGE_REL}}}Relationship")
                      if item.attrib["Id"] == rel_id)
        target = target.lstrip("/")
        worksheet_path = target if target.startswith("xl/") else f"xl/{target}"
        worksheet = ET.fromstring(workbook.read(worksheet_path))
        cells = {cell.attrib["r"]: cell for cell in worksheet.findall(f".//{{{NS_MAIN}}}sheetData/{{{NS_MAIN}}}row/{{{NS_MAIN}}}c")}
        with (ROOT / "data" / "source-data.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        mismatches = []
        for row in rows:
            cell_ref = row["source_cell"]
            cell = cells.get(cell_ref)
            value = cell.find(f"{{{NS_MAIN}}}v") if cell is not None else None
            if value is None or float(value.text) != float(row["relative_affinity_percent"]):
                mismatches.append(cell_ref)
    if mismatches:
        raise SystemExit("Source value mismatch at: " + ", ".join(mismatches))
    print(f"Verified {len(rows)} CSV values against Figure 3d source cells; workbook SHA-256 {digest}")


if __name__ == "__main__":
    main()
