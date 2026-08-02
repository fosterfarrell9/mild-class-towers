"""Sweep-pipeline loader, vendored verbatim for the public tree.

Only `load_values` (with its `LABEL_POINTS` table) is used here; it
parses a six-character matrix table in the sweep's TSV layout.
"""

from __future__ import annotations

import csv
from pathlib import Path

LABEL_POINTS = {
    "x1": (1, 0, 0),
    "x2": (0, 1, 0),
    "x3": (0, 0, 1),
    "x1+x2+x3": (1, 1, 1),
    "x1+x2": (1, 1, 0),
    "x1+x3": (1, 0, 1),
}


def load_values(matrix_file: Path):
    with matrix_file.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 36:
        raise AssertionError(f"{matrix_file}: expected 36 rows, got {len(rows)}")
    primary, doubled, norm_only = {}, {}, {}
    for row in rows:
        target = doubled if row["doubled"] == "1" else primary
        label = row["label"]
        target.setdefault(label, [[0] * 3 for _ in range(3)])
        if row["doubled"] == "0":
            norm_only.setdefault(label, [[0] * 3 for _ in range(3)])
        column = int(row["input"]) - 1
        for output in range(3):
            target[label][output][column] = int(row[f"d{output + 1}"])
            if row["doubled"] == "0":
                norm_only[label][output][column] = int(
                    row[f"norm_only{output + 1}"])
    by_point = {point: primary[label] for label, point in LABEL_POINTS.items()}
    norm_by_point = {point: norm_only[label]
                     for label, point in LABEL_POINTS.items()}
    return rows, primary, doubled, by_point, norm_by_point
