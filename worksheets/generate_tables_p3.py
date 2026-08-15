#!/usr/bin/env python3
"""Generate Part W6 (tables-p3.tex) from the result records.

For every mild field of the range |D_K| < 2^30 the tables list the
discriminant, the class group (from the certified certificate
header), the rank of delta (from the cone-criterion summaries), and
the source of mildness: a transverse cone point with its degree, an
Anick witness, or a terminating Hilbert-series computation.  The
script asserts the published counts (151 + 227 + 481 = 859) and
fails on any inconsistency between the records.

Run from the repository root or from worksheets/:

    python3 worksheets/generate_tables_p3.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
P3 = ROOT / "examples" / "p3"

ANICK = {-211248887, -263780072}
BLOCKS = [
    ("Table 3: the mild fields with $|D_K|<2^{28}$", 0, 2**28, 151),
    ("Table 4: the mild fields with $2^{28}\\le|D_K|<2^{29}$",
     2**28, 2**29, 227),
    ("Table 5: the mild fields with $2^{29}\\le|D_K|<2^{30}$",
     2**29, 2**30, 481),
]


def bucket(disc: int) -> str:
    return f"{abs(disc) // 10**7:03d}"


def load_records(path: Path) -> list[dict]:
    data = json.loads(path.read_text())
    return data if isinstance(data, list) else data["records"]


def class_group(disc: int) -> str:
    certificate = (P3 / "certificates" / bucket(disc)
                   / f"K-{abs(disc)}-p3" / "certificate.gp")
    header = certificate.read_text().splitlines()[0]
    match = re.search(r"\[\[([0-9, ]+)\]", header)
    if not match:
        raise SystemExit(f"no class group in header: {certificate}")
    divisors = match.group(1).replace(" ", "")
    return f"$[{divisors}]$"


def main() -> int:
    cone: dict[int, dict] = {}
    for name in ("cone-criterion/summary.json",
                 "cone-criterion/summary-slice2.json",
                 "cone-criterion/summary-block23-independent.json"):
        for record in load_records(P3 / name):
            cone[record["discriminant"]] = record

    strongly_free: set[int] = set()
    for name in ("results/strong-freeness-block-001.json",
                 "results/slice2-strong-freeness.json",
                 "results/block23-strong-freeness.json"):
        for record in load_records(P3 / name):
            if record["verdict"] == "STRONGLY_FREE":
                strongly_free.add(record["discriminant"])

    rows_by_block: list[list[str]] = []
    total = 0
    for _, lower, upper, expected in BLOCKS:
        mild = []
        for disc, record in cone.items():
            if not lower <= abs(disc) < upper:
                continue
            degree = record["min_transverse_cone_degree"]
            if degree is not None:
                mild.append((disc, f"transverse point, degree ${degree}$"))
            elif disc in ANICK:
                mild.append((disc, "Anick witness"))
            elif disc in strongly_free:
                mild.append((disc, "Hilbert series"))
        if len(mild) != expected:
            raise SystemExit(
                f"count mismatch in [{lower},{upper}): "
                f"{len(mild)} != {expected}")
        mild.sort(key=lambda pair: abs(pair[0]))
        rows = []
        for index, (disc, source) in enumerate(mild, start=total + 1):
            rank = cone[disc]["bockstein_rank"]
            rows.append(f"${index}$ & ${disc}$ & {class_group(disc)} & "
                        f"${rank}$ & {source}\\\\")
        rows_by_block.append(rows)
        total += len(mild)
    assert total == 859, total

    intro = r"""The $859$ mild fields at $p=3$ with $|D_K|<2^{30}$, with the class
group, the rank of $\delta$, and the source of mildness; the three
tables follow the blocks of the computation.  Additional records:
the Anick witnesses in
\nolinkurl{examples/p3/results/block-witnesses-002.json}; the
Gr\"obner verdicts, including the $26$ fields that are provably not
strongly free and the single field with a rank-two tensor, in
\nolinkurl{examples/p3/results/strong-freeness-block-001.json},
\nolinkurl{examples/p3/results/slice2-strong-freeness.json}, and
\nolinkurl{examples/p3/results/block23-strong-freeness.json} ---
every decided verdict was produced by both engines (appendix~C.3
of the paper); the cone-criterion
reports with the transverse cone point of every criterion field in
\nolinkurl{examples/p3/cone-criterion/}.  The drivers are
\nolinkurl{examples/p3/strong-freeness/block_sweep.py},
\nolinkurl{tools/strong_freeness_singular.py}, and
\nolinkurl{tools/p3_cone_criterion.py}.  The exhaustive Anick
searches over the two subfamilies of Section~5.3 of the paper ran
over all $11\,232$ changes of variables over $\Fthree$ and all
$663\,390$ ordered projective bases over $\Fnine$."""

    table_head = (
        "\\begin{center}\\small\n"
        "\\setlength{\\tabcolsep}{5pt}%\n"
        "\\begin{longtable}{r r l c l}\n"
        "\\toprule\n"
        "& $D_K$ & $\\Cl(K)$ & $\\rk\\delta$ & mild via\\\\\n"
        "\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n")

    parts = ["\\section*{Part W6: the $p=3$ tables}", "", intro, ""]
    for (title, *_), rows in zip(BLOCKS, rows_by_block):
        parts += ["\\subsection*{" + title + "}", "", table_head]
        parts += rows
        parts += ["\\end{longtable}", "\\end{center}", ""]
    body = "\n".join(parts)

    preamble = re.split(r"\\begin\{document\}",
                        (HERE / "worksheets.tex").read_text())[0]
    document = (
        preamble + "\\begin{document}\n\n\\begin{center}\n"
        "{\\Large Worksheet Part W6: the $p=3$ tables}\\\\[2mm]\n"
        "accompanying \\emph{Mild $p$-class tower groups of imaginary "
        "quadratic\nfields}\n\\end{center}\n\n"
        "\\medskip\n\\noindent\n"
        "This document is generated from the result records by\n"
        "\\nolinkurl{generate_tables_p3.py}; references of the form\n"
        "Section~5.3 refer to the paper.\n\n"
        + body + "\n\\end{document}\n")
    (HERE / "tables-p3.tex").write_text(document)
    print(f"tables-p3.tex written: {total} fields, "
          f"{sum(len(r) for r in rows_by_block)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
