#!/usr/bin/env python3
"""Render the replicate throughput chart from the NVIDIA OpenMM run record."""

from __future__ import annotations

import json
import math
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "data" / "nvidia-engine-preflight.json"
OUTPUT = ROOT / "assets" / "nvidia-t4-preflight.svg"
WIDTH = 1000
HEIGHT = 360
LEFT = 210
RIGHT = 835
TOP = 142
ROW_GAP = 48
BAR_HEIGHT = 16
PURPLE = "#6c53a2"
PLUM = "#302541"
MUTED = "#77737e"
MINT = "#cde6aa"


def main() -> None:
    record = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = record["measurements"]
    summary = record["summary"]
    maximum = max(row["ns_per_day"] for row in rows)
    scale = math.ceil(maximum / 50.0) * 50.0
    plot_width = RIGHT - LEFT
    run_id = escape(record["run"]["id"])
    device = escape(record["protocol"]["device_name"])
    atom_count = f"{record['system']['atom_count']:,}"
    median = summary["median_ns_per_day"]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-labelledby="title desc">',
        '<title id="title">OpenMM OpenCL benchmark throughput on an NVIDIA GPU</title>',
        f'<desc id="desc">Three 100-picosecond timed runs on {device} using a generic solvated 1UBQ system with {atom_count} atoms. Median throughput {median:.1f} nanoseconds per day. This is not an NNC2215 or PaG simulation.</desc>',
        f'<rect width="{WIDTH}" height="{HEIGHT}" rx="14" fill="#ffffff"/>',
        f'<text x="32" y="38" fill="{PURPLE}" font-family="Arial,sans-serif" font-size="11" font-weight="700" letter-spacing="1.6">NVIDIA GPU · ENGINE PREFLIGHT</text>',
        f'<text x="32" y="72" fill="{PLUM}" font-family="Arial,sans-serif" font-size="23" font-weight="700">OpenMM OpenCL throughput</text>',
        f'<text x="32" y="99" fill="{MUTED}" font-family="Arial,sans-serif" font-size="12">{device} · {atom_count} atoms · 3 independent 100 ps windows</text>',
        f'<text x="965" y="47" text-anchor="end" fill="{MUTED}" font-family="Arial,sans-serif" font-size="10">{run_id}</text>',
        f'<rect x="805" y="62" width="160" height="54" rx="10" fill="#f2f0f7"/>',
        f'<text x="822" y="84" fill="{MUTED}" font-family="Arial,sans-serif" font-size="9" font-weight="700" letter-spacing="1">MEDIAN</text>',
        f'<text x="822" y="105" fill="{PLUM}" font-family="Arial,sans-serif" font-size="18" font-weight="700">{median:.1f} ns/day</text>',
    ]

    for tick in range(0, int(scale) + 1, 50):
        x = LEFT + plot_width * tick / scale
        parts.append(
            f'<line x1="{x:.1f}" y1="{TOP - 16}" x2="{x:.1f}" y2="{TOP + ROW_GAP * 2 + 20}" stroke="#e6e3e9" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{TOP + ROW_GAP * 2 + 40}" text-anchor="middle" fill="{MUTED}" font-family="Arial,sans-serif" font-size="10">{tick}</text>'
        )

    for index, row in enumerate(rows):
        y = TOP + index * ROW_GAP
        width = plot_width * row["ns_per_day"] / scale
        label = f"Run {index + 1} · seed {row['seed']}"
        parts.append(
            f'<text x="32" y="{y + 4}" fill="{PLUM}" font-family="Arial,sans-serif" font-size="12" font-weight="700">{escape(label)}</text>'
        )
        parts.append(
            f'<rect x="{LEFT}" y="{y - BAR_HEIGHT / 2}" width="{plot_width}" height="{BAR_HEIGHT}" rx="8" fill="#f0eef3"/>'
        )
        parts.append(
            f'<rect x="{LEFT}" y="{y - BAR_HEIGHT / 2}" width="{width:.1f}" height="{BAR_HEIGHT}" rx="8" fill="{PURPLE}"/>'
        )
        parts.append(
            f'<circle cx="{LEFT + width:.1f}" cy="{y:.1f}" r="5" fill="{MINT}" stroke="{PLUM}" stroke-width="1.5"/>'
        )
        parts.append(
            f'<text x="{min(LEFT + width + 14, 955):.1f}" y="{y + 4}" fill="{PLUM}" font-family="Arial,sans-serif" font-size="12" font-weight="700">{row["ns_per_day"]:.1f}</text>'
        )

    parts.extend(
        [
            f'<text x="{LEFT + plot_width / 2:.1f}" y="{HEIGHT - 17}" text-anchor="middle" fill="{MUTED}" font-family="Arial,sans-serif" font-size="10">SIMULATED THROUGHPUT · NANOSECONDS PER DAY</text>',
            "</svg>",
        ]
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)} for {run_id}")


if __name__ == "__main__":
    main()
