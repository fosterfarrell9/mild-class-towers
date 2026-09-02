#!/usr/bin/env python3
"""Generate W5-tables-p5.tex, the unified p = 5 table (Part W5).

One table for all 204 fields of the census, with the same columns for
every field: class group, number of F_5-rational points of Sigma_D at
which rk D_x = 1, number of these which are transverse, and the least
degree of an extension of F_5 carrying a transverse point.  The counts
are read from records/p5/transverse-rank-one/certificates.json and
cross-checked against records/p5/transverse-rank-one/report.txt; the
class groups are read from the arithmetic certificates
certificates/p5/K-*-p5/certificate.gp.  The field D_K = -781922404,
whose scheme carries no transverse point over any extension, is listed
with a footnote pointing to Part W7.

Run from the repository root:  python3 worksheets/generate_tables_p5.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
JSON = ROOT / "records" / "p5" / "transverse-rank-one" / "certificates.json"
REPORT = ROOT / "records" / "p5" / "transverse-rank-one" / "report.txt"
CERTS = ROOT / "certificates" / "p5"
OUT = ROOT / "worksheets" / "W5-tables-p5.tex"

EXCEPTION = -781922404


def radicand(disc: int) -> int:
    """Directory number of the field: |D| if D is odd, else |D|/4."""
    a = abs(disc)
    return a if disc % 2 != 0 else a // 4


def class_group(disc: int) -> list[int]:
    path = CERTS / f"K-{radicand(disc)}-p5" / "certificate.gp"
    head = path.read_text().split("\n", 1)[0]
    if str(disc) not in head:
        raise SystemExit(f"{path}: discriminant {disc} not in header")
    m = re.search(r"\[\[(\d+(?:,\s*\d+)*)\]", head)
    if not m:
        raise SystemExit(f"{path}: no class group in header")
    return [int(x) for x in m.group(1).split(",")]


def load_fields() -> list[dict]:
    fields = []
    for e in json.load(open(JSON)):
        pts = e["rational_points"]
        if any(p["rank"] != 1 for p in pts):
            raise SystemExit(f"{e['disc']}: rational point of rank != 1")
        fields.append({
            "disc": e["disc"],
            "rk1": len(pts),
            "trans": sum(1 for p in pts if p["transverse"]),
            "mindeg": e["minimal_degree"],
            "cl": class_group(e["disc"]),
        })
    fields.sort(key=lambda f: abs(f["disc"]))
    return fields


def crosscheck_report(fields: list[dict]) -> None:
    rows = {}
    for line in open(REPORT):
        m = re.match(r"(-\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(YES|NO)", line)
        if m:
            rows[int(m.group(1))] = (int(m.group(2)), int(m.group(3)),
                                     None if m.group(4) in {"-", "None"}
                                     else int(m.group(4)))
    if len(rows) != len(fields):
        raise SystemExit(f"report.txt: {len(rows)} rows, expected {len(fields)}")
    for f in fields:
        got = (f["rk1"], f["trans"], f["mindeg"])
        if rows[f["disc"]] != got:
            raise SystemExit(f"{f['disc']}: json {got} != report {rows[f['disc']]}")


HEADER = r"""\documentclass[11pt]{article}
\usepackage[margin=2.6cm]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{booktabs,longtable}
\newcommand{\Ffive}{\mathbb F_5}
\newcommand{\Cl}{\operatorname{Cl}}
\newcommand{\rk}{\operatorname{rk}}
\usepackage{url}
\providecommand{\nolinkurl}[1]{\url{#1}}
\begin{document}

\begin{center}
{\Large Worksheet Part W5: the $p=5$ table}\\[2mm]
accompanying \emph{Mild $p$-class tower groups of imaginary quadratic
fields}
\end{center}

\medskip
\noindent
References of the form Definition~3.17 or Part~W7 refer to the
statements of the paper and to the parts of this repository.

\section*{Part W5: The $p=5$ table}

The $204$ imaginary quadratic fields with $5$-class rank three and
$|D_K|<2^{30}$, one line per field.  The columns record the class
group, the number of $\Ffive$-rational points of the norm-degeneracy
scheme $\Sigma_D$ (Definition~3.17 of the paper) at which
$\rk D_x=1$, the number of these which are transverse, and the least
degree of an extension of $\Ffive$ over which a transverse point
exists.  The defining polynomials of the extensions used for the
fields without a rational transverse point, the coordinates of all
transverse points, and the corresponding transversality matrices are
recorded for every field in
\nolinkurl{records/p5/transverse-rank-one/}
(\nolinkurl{report.txt} human-readable,
\nolinkurl{certificates.json} machine-readable); the finite-field
calculations underlying the table are re-run deterministically by
\nolinkurl{tools/transverse_rank_one.py}.  For the field
$D_K=-781922404$, marked ${}^{*}$, the scheme $\Sigma_D$ is a single
non-reduced closed point of degree two, so no transverse point exists
over any extension; the field is mild by a word count after a change
of variables over $\mathbb F_{25}$, recorded in Part~W7.

\begin{center}\small
\setlength{\tabcolsep}{4.5pt}%
\begin{longtable}{r r l c c c}
\toprule
& $D_K$ & $\Cl(K)$ & rank-one over $\Ffive$ & transverse over $\Ffive$
& \shortstack{least transverse\\degree} \\
\midrule
\endhead
\bottomrule
\endfoot
"""

FOOTER = r"""\end{longtable}
\end{center}

\end{document}
"""


def main() -> int:
    fields = load_fields()
    crosscheck_report(fields)
    lines = []
    for i, f in enumerate(fields, 1):
        cl = "[" + ",".join(str(x) for x in f["cl"]) + "]"
        deg = str(f["mindeg"]) if f["mindeg"] is not None else "---${}^{*}$"
        lines.append(f"${i}$ & ${f['disc']}$ & ${cl}$ & ${f['rk1']}$ &"
                     f" ${f['trans']}$ & {deg if f['mindeg'] is None else f'${deg}$'}\\\\")
    OUT.write_text(HEADER + "\n".join(lines) + "\n" + FOOTER)
    n_trans = sum(1 for f in fields if f["mindeg"] is not None)
    n_rat = sum(1 for f in fields if f["mindeg"] == 1)
    print(f"{len(fields)} Zeilen geschrieben, {n_trans} mit transversalem"
          f" Punkt ({n_rat} rational), 1 Ausnahme")
    return 0


if __name__ == "__main__":
    sys.exit(main())
