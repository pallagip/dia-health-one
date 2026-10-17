#!/usr/bin/env python3
"""Audit the published PaG receptor-assay text and source-data availability."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parent
XML_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13270365/fullTextXML"
ARTICLE_URL = "https://pubs.acs.org/doi/10.1021/acsptsci.6c00069"
PMC_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13270365/"
SI_AUDIT_PATH = ROOT / "data" / "pag-si-source-audit.json"
OUTPUT_PATH = ROOT / "data" / "pag-assay-source-audit.json"


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def element_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def fetch_article_xml() -> bytes:
    request = Request(XML_URL, headers={"User-Agent": "DiaHealthOneResearchAudit/1.0"})
    with urlopen(request, timeout=30) as response:
        return response.read()


def find_text(root: ET.Element, tag: str, predicate) -> str:
    for element in root.iter():
        if local_name(element.tag) == tag:
            text = element_text(element)
            if predicate(text):
                return text
    raise ValueError(f"Could not find the expected {tag} source passage")


def main() -> None:
    xml_bytes = fetch_article_xml()
    root = ET.fromstring(xml_bytes)
    ids = {
        element.attrib.get("pub-id-type"): element_text(element)
        for element in root.iter()
        if local_name(element.tag) == "article-id"
    }
    if ids.get("doi") != "10.1021/acsptsci.6c00069":
        raise ValueError(f"Unexpected article DOI in source XML: {ids.get('doi')!r}")

    title = next(
        element_text(element)
        for element in root.iter()
        if local_name(element.tag) == "article-title"
    )
    license_text = find_text(root, "license", lambda text: "CC-BY-NC-ND 4.0" in text)
    assay_text = find_text(
        root,
        "p",
        lambda text: "0.9%" in text and "5.4%" in text and "70 μmol/L" in text,
    )
    match = re.search(r"from\s+(\d+(?:\.\d+)?)%\s+to\s+(\d+(?:\.\d+)?)%", assay_text)
    if match is None:
        raise ValueError("Could not parse the two article-text receptor-affinity values")
    no_pag_affinity = float(match.group(1))
    pag_affinity = float(match.group(2))

    figure_caption = find_text(
        root,
        "fig",
        lambda text: "Glucose-sensitive binding of NNC2215" in text and "hIR-A" in text,
    )
    if "mean and SD" not in figure_caption or "N = 3" not in figure_caption:
        raise ValueError("Figure 3 caption did not match the reviewed mean/SD, N=3 description")

    p_value_text = find_text(root, "p", lambda text: "p < 0.0001" in text and "two-way" in text)
    p_value = re.search(r"p\s*<\s*0\.0001", p_value_text)
    if p_value is None:
        raise ValueError("Could not confirm the article's reported two-way ANOVA value")

    si_bytes = SI_AUDIT_PATH.read_bytes()
    si_audit = json.loads(si_bytes)
    if si_audit["public_record_availability"]["figshare_item_file_count"] != 1:
        raise ValueError("Supporting Information inventory changed; review before publishing")

    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    script_bytes = Path(__file__).read_bytes()
    ratio = pag_affinity / no_pag_affinity
    record = {
        "schema": "dia-health-one.pag-assay-source-audit.v1",
        "recorded_at_utc": generated_at,
        "run": {
            "script": "research/audit_pag_assay_source.py",
            "script_sha256": hashlib.sha256(script_bytes).hexdigest(),
            "article_full_text_xml_sha256": hashlib.sha256(xml_bytes).hexdigest(),
            "supporting_information_audit": "research/data/pag-si-source-audit.json",
            "supporting_information_audit_sha256": hashlib.sha256(si_bytes).hexdigest(),
        },
        "source": {
            "title": title,
            "doi": ids["doi"],
            "pmcid": ids["pmcid"],
            "article_url": ARTICLE_URL,
            "full_text_xml_url": XML_URL,
            "pmc_url": PMC_URL,
            "license": {
                "name": "CC BY-NC-ND 4.0",
                "url": "https://creativecommons.org/licenses/by-nc-nd/4.0/",
                "source_xml_statement_verified": "CC-BY-NC-ND 4.0" in license_text,
            },
        },
        "assay": {
            "target": "NNC2215 binding to human insulin receptor isoform A (hIR-A)",
            "matrix": "1.5% human serum albumin",
            "pag_condition_umol_l": 70,
            "figure_caption_summary": "Figure 3 reports mean and SD error bars from three independent experiments.",
            "reported_two_way_anova_p": "<0.0001",
            "p_value_recalculated": False,
        },
        "exact_article_text_values": [
            {
                "glucose_mmol_l": 0,
                "pag_umol_l": 0,
                "relative_affinity_percent_of_human_insulin": no_pag_affinity,
                "evidence_status": "value quoted in article text",
            },
            {
                "glucose_mmol_l": 0,
                "pag_umol_l": 70,
                "relative_affinity_percent_of_human_insulin": pag_affinity,
                "evidence_status": "value quoted in article text",
            },
        ],
        "derived_from_exact_text_values": {
            "pag_to_no_pag_fold_change_at_zero_glucose": ratio,
            "absolute_change_percentage_points": pag_affinity - no_pag_affinity,
            "source_uncertainty_recalculated": False,
        },
        "source_data_availability": {
            "figure_3_has_visual_mean_sd_points": True,
            "machine_readable_curve_means_found_in_reviewed_article_or_si": False,
            "machine_readable_numeric_sd_values_found_in_reviewed_article_or_si": False,
            "individual_replicate_values_found_in_reviewed_article_or_si": False,
            "supporting_information_file_count": si_audit["public_record_availability"]["figshare_item_file_count"],
            "supporting_information_file": si_audit["source"]["file"]["name"],
            "supporting_information_pages": si_audit["source"]["file"]["page_count"],
        },
        "public_demo_reuse_decision": {
            "figure_embedded_or_adapted": False,
            "figure_digitized_for_public_dashboard": False,
            "reason": "The reviewed source is licensed CC BY-NC-ND 4.0 and does not provide a machine-readable point table. The demo links the publisher source and displays only the exact article-text values recorded above.",
        },
        "limits": [
            "The 0.9% and 5.4% values are rounded published means; the underlying replicate values and exact numerical SDs were not found in the reviewed public files.",
            "The reported p-value is the paper's two-way ANOVA result and was not recalculated.",
            "The derived 6.0-fold ratio uses only the two rounded article-text values; it is not a molecular simulation or patient-specific effect.",
        ],
        "clinical_use": False,
    }
    OUTPUT_PATH.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT_PATH), "recorded_at_utc": generated_at, "values_percent": [no_pag_affinity, pag_affinity], "ratio": ratio, "license": record["source"]["license"]["name"], "source_xml_sha256": record["run"]["article_full_text_xml_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
