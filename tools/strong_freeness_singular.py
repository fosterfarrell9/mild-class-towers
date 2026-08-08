#!/usr/bin/env python3
"""Singular/Letterplace driver for the strong-freeness Groebner test.

This is an independent engine for the computation performed by
strong_freeness_gb.py: the truncated two-sided Groebner basis of the cubic
relation ideal is computed by Singular's Letterplace subsystem
(La Scala--Levandovskyy), and only the leading words are passed back; the
normal-word counting and the comparison with the strongly free target
series 1/(1-3z+3z^3) reuse the audited routines of strong_freeness_gb.py.

The characteristic is taken from --prime and enters both the generated
Singular ring and the reduction of the input coefficients.  Note that a
caller importing this module gets its own instance of
strong_freeness_gb; switch the characteristic through ``sing.gb.P`` (or
--prime), never through a separately imported copy, or the ring is built
over the wrong field while the coefficients still look plausible.

The --order option uses the same ascending convention as
strong_freeness_gb.py: --order xyz means x < y < z in the
degree-lexicographic order, so words rich in the last letter lead.  The
generated Singular ring lists the variables in descending order
accordingly.

Termination is detected by the diamond-lemma bound: if every basis element
has degree at most m and the degree bound is at least 2m-1, all overlap
ambiguities have been processed and the basis is a complete Groebner
basis.  The verdicts follow strong_freeness_gb.py:

  STRONGLY_FREE               complete basis, series equals the target
                              (rigorous by the linear-recurrence bound);
  NOT_STRONGLY_FREE           a verified coefficient deviates;
  INCONCLUSIVE_SERIES_MATCHES truncated basis, series verified through
                              the degree bound.

Use --emit-script to print the generated Singular input instead of
running it, for archival or manual runs.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "strong_freeness_gb", HERE / "strong_freeness_gb.py")
assert SPEC is not None and SPEC.loader is not None
gb = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = gb
SPEC.loader.exec_module(gb)


def polynomial_text(row):
    terms = []
    for column, coeff in enumerate(row):
        if coeff % gb.P:
            i, j, k = column // 9, (column // 3) % 3, column % 3
            word = "*".join(gb.LETTERS[t] for t in (i, j, k))
            terms.append(f"{coeff % gb.P}*{word}")
    return " + ".join(terms)


def build_script(matrix, order, bound):
    descending = ",".join(reversed(order))
    return f"""LIB "freegb.lib";
ring r = {gb.P},({descending}),Dp;
def R = freeAlgebra(r, {bound});
setring R;
ideal I = {polynomial_text(matrix[0])},
  {polynomial_text(matrix[1])},
  {polynomial_text(matrix[2])};
option(redSB); option(redTail);
ideal J = twostd(I);
print("GB_LEADS");
for (int q = 1; q <= size(J); q++) {{ print(lead(J[q])); }}
print("GB_END");
exit;
"""


def parse_leads(output: str):
    lines = output.splitlines()
    try:
        start = lines.index("GB_LEADS") + 1
        end = lines.index("GB_END")
    except ValueError:
        raise SystemExit(
            "Singular output did not contain the GB_LEADS/GB_END markers:\n"
            + output[-2000:])
    leads = []
    for line in lines[start:end]:
        word = re.sub(r"[^xyz]", "", line)
        if word:
            leads.append(word)
    return leads


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--result", type=Path)
    source.add_argument("--matrix")
    parser.add_argument("--order", default="xyz")
    parser.add_argument("--prime", type=int, default=gb.P,
                        help="characteristic of the coefficient field; the "
                             "p=3 block drivers pass 3")
    parser.add_argument("--degree-bound", type=int, default=13)
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--emit-script", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    if sorted(args.order) != ["x", "y", "z"]:
        raise SystemExit("--order must be a permutation of xyz")
    gb.P = args.prime
    matrix = (gb.matrix_from_result(args.result) if args.result
              else gb.parse_matrix_text(args.matrix))
    if not gb.rank3(matrix):
        raise SystemExit(f"matrix does not have rank 3 over F_{gb.P}")

    script = build_script(matrix, args.order, args.degree_bound)
    if args.emit_script:
        print(script)
        return 0
    if shutil.which(args.singular) is None:
        raise SystemExit(f"Singular executable not found: {args.singular}")

    completed = subprocess.run(
        [args.singular, "-q"], input=script, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    leads = parse_leads(completed.stdout)
    largest = max(len(w) for w in leads)
    terminated = args.degree_bound >= 2 * largest - 1
    by_degree = {}
    for lead in leads:
        by_degree[len(lead)] = by_degree.get(len(lead), 0) + 1
    print(f"order = deg-lex {args.order} (Singular Dp); "
          f"degree bound = {args.degree_bound}; terminated = {terminated}")
    print("basis leading words by degree: "
          + ", ".join(f"{d}:{c}" for d, c in sorted(by_degree.items())))

    if terminated:
        probe = 1
        series, states = gb.normal_word_counts(leads, probe)
        probe = 2 * (states + 3) + 3
        series, states = gb.normal_word_counts(leads, probe)
        target = gb.target_series(probe)
        verdict = ("STRONGLY_FREE" if series == target
                   else "NOT_STRONGLY_FREE")
        print(f"automaton states = {states}; comparison through degree "
              f"{probe} (rigorous bound 2(s+3))")
    else:
        series, states = gb.normal_word_counts(leads, args.degree_bound)
        target = gb.target_series(args.degree_bound)
        verdict = ("INCONCLUSIVE_SERIES_MATCHES" if series == target
                   else "NOT_STRONGLY_FREE")
        print(f"series verified through degree {args.degree_bound}: "
              f"{'matches target' if series == target else 'DEVIATES'}")
    if series != target:
        first_bad = next(
            n for n in range(len(series)) if series[n] != target[n])
        print(f"first deviation at degree {first_bad}: "
              f"a_n = {series[first_bad]}, target {target[first_bad]}")
    print(f"VERDICT: {verdict}")
    if args.json:
        print(json.dumps({
            "engine": "singular-letterplace",
            "order": args.order,
            "degree_bound": args.degree_bound,
            "terminated": terminated,
            "leading_words": sorted(leads),
            "verdict": verdict,
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
