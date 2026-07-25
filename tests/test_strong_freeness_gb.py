#!/usr/bin/env python3
"""Tests for the noncommutative Groebner strong-freeness tool."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "strong_freeness_gb.py"
SPEC = importlib.util.spec_from_file_location("strong_freeness_gb", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
gb = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = gb
SPEC.loader.exec_module(gb)


def run(matrix, order="xyz", max_degree=10):
    key = gb.make_order_key(order)
    relations = gb.rows_to_polynomials(matrix)
    basis, processed, terminated = gb.complete(relations, key, max_degree)
    leads = sorted(basis, key=key)
    return basis, leads, processed, terminated


# 1. The transformed relations of the paper's principal example have the
#    combinatorially free leading words zzx, zyy, zyx; the completion must
#    terminate immediately and certify strong freeness.
UT_PHI = gb.parse_matrix_text(
    "[0,2,0,1,4,0,0,2,1,2,2,3,4,0,0,2,0,0,0,3,3,0,0,0,1,0,0;"
    "0,3,4,4,2,0,2,1,0,3,1,4,2,0,1,1,3,0,4,4,0,0,1,0,0,0,0;"
    "0,1,2,3,2,1,1,4,0,1,1,0,2,0,0,4,0,0,2,0,0,1,0,0,0,0,0]")
basis, leads, processed, terminated = run(UT_PHI)
assert terminated, "completion must terminate for combinatorially free leads"
assert leads == ["zyx", "zyy", "zzx"], leads
series, states = gb.normal_word_counts(leads, 2 * (4 + 3) + 3)
assert states == 4, states
assert series == gb.target_series(len(series) - 1), "series must match target"

# 2. The monomial relations xxx, xxy, xxz form their own Groebner basis but
#    are not strongly free: the quotient has 66 > 63 words in degree 4.
MONOMIAL = [[0] * 27 for _ in range(3)]
for column in range(3):
    MONOMIAL[column][column] = 1  # columns 1..3 are xxx, xxy, xxz
basis, leads, processed, terminated = run(MONOMIAL)
assert terminated
assert leads == ["xxx", "xxy", "xxz"], leads
series, states = gb.normal_word_counts(leads, 6)
assert series[:5] == [1, 3, 9, 24, 66], series[:5]
assert gb.target_series(4)[4] == 63

# 3. Reduction sanity: reducing a relation by itself gives zero.
key = gb.make_order_key("xyz")
relations = gb.rows_to_polynomials(UT_PHI)
lead = gb.leading_word(relations[0], key)
reduced = gb.reduce_poly(relations[0], {lead: gb.monic(relations[0], key)}, key)
assert reduced == {}, reduced

# 4. Automaton fallback: single leading word zzz leaves 3^n - #(words
#    containing zzz); check degree 3 and 4 counts directly.
series, states = gb.normal_word_counts(["zzz"], 4)
assert series == [1, 3, 9, 26, 76], series

print("STRONG_FREENESS_GB_TEST PASS")
