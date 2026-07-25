#!/usr/bin/env python3
"""Noncommutative Groebner-basis test for strong freeness of cubic relations.

Input is a 3 x 27 matrix T over F_5 whose rows are the coordinate vectors of
three cubic relations in F_5<x,y,z> (word order: X_i X_j X_k with k fastest,
as in the paper and in the result.gp records).  The rows must span a
three-dimensional relation space R_3.

The script fixes a degree-lexicographic monomial order (letter order
configurable), runs Bergman's diamond-lemma completion degree by degree, and
decides:

  * If the completion TERMINATES (all overlap ambiguities up to composite
    degree 2m-1 resolve, where m is the largest degree of a basis element),
    the set of leading words is a finite Groebner leading set.  The Hilbert
    series of A/(R_3) then equals the growth series of the words avoiding
    the leading set, which is computed exactly through an Aho--Corasick
    automaton, and comparison with the strongly free target series

        1/(1 - 3z + 3z^3)

    becomes a finite, rigorous check: both sequences satisfy linear
    recurrences of order at most s+3 (s = number of automaton states), so
    agreement in degrees <= 2(s+3) implies agreement in every degree.
    The verdict is then STRONGLY_FREE or NOT_STRONGLY_FREE.

  * If the completion does not terminate before --max-degree, the verdict is
    INCONCLUSIVE, but the Hilbert series of A/(R_3) is still rigorously
    verified against the target in every completed degree (the leading words
    of the partial basis determine dim (A/(R_3))_n exactly for all n up to
    the last fully processed degree).

Failure to certify strong freeness proves nothing about mildness; a verdict
NOT_STRONGLY_FREE proves that the group is not mild with respect to the
Zassenhaus filtration; STRONGLY_FREE proves that it is mild.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

P = 5
LETTERS = "xyz"


# ----------------------------------------------------------------------
# polynomials: dict mapping words (strings over "xyz") to nonzero coeffs


def poly_normalize(poly: dict) -> dict:
    return {w: c % P for w, c in poly.items() if c % P}


def poly_sub_scaled(target: dict, factor: int, source: dict) -> None:
    """target -= factor * source, in place."""
    for word, coeff in source.items():
        value = (target.get(word, 0) - factor * coeff) % P
        if value:
            target[word] = value
        else:
            target.pop(word, None)


def leading_word(poly: dict, key) -> str:
    return max(poly, key=key)


# ----------------------------------------------------------------------
# monomial order: degree first, then lexicographic in a letter permutation


def make_order_key(letter_order: str):
    rank = {ch: i for i, ch in enumerate(letter_order)}

    def key(word: str):
        return (len(word), [rank[ch] for ch in word])

    return key


# ----------------------------------------------------------------------
# reduction and completion


def find_division(word: str, leads: dict):
    """Return (position, lead) for the first occurrence of a basis leading
    word inside `word`, or None."""
    for lead in leads:
        position = word.find(lead)
        if position >= 0:
            return position, lead
    return None


def reduce_poly(poly: dict, basis: dict, key) -> dict:
    """Fully reduce poly modulo the basis (leading word -> monic poly).

    Words already verified irreducible stay irreducible while the basis is
    fixed, so each word is tested at most once; a lazy max-heap orders the
    candidates."""
    import heapq

    def heap_entry(word):
        k = key(word)
        return (-k[0], [-r for r in k[1]], word)

    poly = dict(poly)
    heap = [heap_entry(w) for w in poly]
    heapq.heapify(heap)
    normal = set()
    while heap:
        word = heapq.heappop(heap)[2]
        if word not in poly or word in normal:
            continue
        hit = find_division(word, basis)
        if hit is None:
            normal.add(word)
            continue
        position, lead = hit
        left, right = word[:position], word[position + len(lead):]
        factor = poly[word]
        shifted = {left + w + right: c for w, c in basis[lead].items()}
        poly_sub_scaled(poly, factor, shifted)
        for fresh in shifted:
            if fresh in poly and fresh not in normal:
                heapq.heappush(heap, heap_entry(fresh))
    return poly


def monic(poly: dict, key) -> dict:
    lead = leading_word(poly, key)
    inverse = pow(poly[lead], P - 2, P)
    return {w: (c * inverse) % P for w, c in poly.items()}


def overlap_ambiguities(w1: str, w2: str):
    """Yield composite words for proper overlaps: suffix of w1 = prefix of
    w2 of length k, 1 <= k < min(|w1|, |w2|)."""
    for k in range(1, min(len(w1), len(w2))):
        if w1[len(w1) - k:] == w2[:k]:
            yield w1 + w2[k:], k


def complete(relations, key, max_degree, progress=None):
    """Degree-truncated diamond-lemma completion.

    Returns (basis, processed_degree, terminated) where basis maps leading
    words to monic fully reduced polynomials, and processed_degree is the
    largest degree D such that every ambiguity of composite degree <= D has
    been resolved.  If `progress` is given it is called as
    progress(degree, basis) after each completed degree.
    """
    basis = {}
    pending = []  # homogeneous polynomials awaiting insertion
    for relation in relations:
        pending.append(poly_normalize(relation))

    def insert_all():
        changed = True
        while changed:
            changed = False
            nonlocal pending
            queue, pending = pending, []
            for poly in queue:
                reduced = reduce_poly(poly, basis, key)
                if not reduced:
                    continue
                reduced = monic(reduced, key)
                lead = leading_word(reduced, key)
                # inter-reduce: retire basis elements whose lead contains
                # the new lead as a subword
                retired = [
                    old for old in basis
                    if old != lead and lead in old
                ]
                basis[lead] = reduced
                for old in retired:
                    element = basis.pop(old)
                    pending.append(element)
                changed = True

    insert_all()
    degree = 3
    while degree < max_degree:
        degree += 1
        new_elements = []
        leads = sorted(basis, key=key)
        for w1 in leads:
            for w2 in leads:
                for composite, _ in overlap_ambiguities(w1, w2):
                    if len(composite) != degree:
                        continue
                    g1, g2 = basis.get(w1), basis.get(w2)
                    if g1 is None or g2 is None:
                        continue
                    right = composite[len(w1):]
                    left = composite[:len(composite) - len(w2)]
                    s_poly = {w + right: c for w, c in g1.items()}
                    poly_sub_scaled(
                        s_poly, 1,
                        {left + w: c for w, c in g2.items()})
                    remainder = reduce_poly(s_poly, basis, key)
                    if remainder:
                        new_elements.append(remainder)
        if new_elements:
            pending.extend(new_elements)
            insert_all()
        if progress is not None:
            progress(degree, basis)
        largest = max(len(w) for w in basis)
        if degree >= 2 * largest - 1:
            return basis, degree, True
    return basis, max_degree, False


# ----------------------------------------------------------------------
# normal-word counting through an Aho--Corasick style automaton


def build_automaton(leads):
    """States are proper prefixes of leading words that contain no leading
    word; transitions follow the longest-suffix rule."""
    prefixes = {""}
    for lead in leads:
        for i in range(1, len(lead)):
            prefixes.add(lead[:i])
    states = sorted(
        p for p in prefixes
        if find_division(p, dict.fromkeys(leads)) is None)
    index = {s: i for i, s in enumerate(states)}
    transitions = []
    for state in states:
        row = []
        for ch in LETTERS:
            word = state + ch
            if any(word.endswith(lead) for lead in leads):
                row.append(None)  # forbidden
            else:
                for start in range(len(word) + 1):
                    if word[start:] in index:
                        row.append(index[word[start:]])
                        break
        transitions.append(row)
    return states, transitions


def normal_word_counts(leads, length):
    states, transitions = build_automaton(leads)
    counts = [0] * len(states)
    counts[states.index("")] = 1
    series = [1]
    for _ in range(length):
        fresh = [0] * len(states)
        for state, value in enumerate(counts):
            if not value:
                continue
            for target in transitions[state]:
                if target is not None:
                    fresh[target] += value
        counts = fresh
        series.append(sum(counts))
    return series, len(states)


def target_series(length):
    a = [1, 3, 9]
    for n in range(3, length + 1):
        a.append(3 * a[n - 1] - 3 * a[n - 3])
    return a[:length + 1]


# ----------------------------------------------------------------------
# input handling


def parse_matrix_text(text: str):
    rows = text.strip().lstrip("[").rstrip("]").split(";")
    matrix = [[int(v) % P for v in row.split(",")] for row in rows]
    if len(matrix) != 3 or any(len(r) != 27 for r in matrix):
        raise SystemExit("matrix must be 3 x 27")
    return matrix


def matrix_from_result(path: Path):
    text = path.read_text()
    match = re.search(
        r'"cubic_relation_matrix",\s*\[(.*?)\]\]', text, re.S)
    if not match:
        raise SystemExit(f"{path}: no cubic_relation_matrix entry")
    return parse_matrix_text("[" + match.group(1) + "]")


def rows_to_polynomials(matrix):
    words = [
        LETTERS[i] + LETTERS[j] + LETTERS[k]
        for i in range(3) for j in range(3) for k in range(3)
    ]
    return [
        {words[column]: value for column, value in enumerate(row) if value}
        for row in matrix
    ]


def rank3(matrix):
    m = [row[:] for row in matrix]
    rank = 0
    for col in range(27):
        pivot = next(
            (r for r in range(rank, 3) if m[r][col] % P), None)
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        inv = pow(m[rank][col], P - 2, P)
        m[rank] = [(v * inv) % P for v in m[rank]]
        for r in range(3):
            if r != rank and m[r][col] % P:
                f = m[r][col]
                m[r] = [(m[r][c] - f * m[rank][c]) % P for c in range(27)]
        rank += 1
    return rank == 3


# ----------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--result", type=Path,
                        help="result.gp file with a cubic_relation_matrix")
    source.add_argument("--matrix",
                        help="matrix literal [r1;r2;r3], rows of length 27")
    parser.add_argument("--order", default="xyz",
                        help="letter order for degree-lex (default xyz)")
    parser.add_argument("--max-degree", type=int, default=14)
    parser.add_argument("--json", action="store_true",
                        help="emit a machine-readable summary line")
    args = parser.parse_args(argv)

    if sorted(args.order) != ["x", "y", "z"]:
        raise SystemExit("--order must be a permutation of xyz")
    matrix = (matrix_from_result(args.result) if args.result
              else parse_matrix_text(args.matrix))
    if not rank3(matrix):
        raise SystemExit("matrix does not have rank 3 over F_5")

    key = make_order_key(args.order)
    relations = rows_to_polynomials(matrix)

    def progress(degree, current):
        leads_now = sorted(current, key=key)
        series, _ = normal_word_counts(leads_now, degree)
        match = series == target_series(degree)
        print(f"  degree {degree}: basis size {len(current)}; "
              f"series {'matches' if match else 'DEVIATES'} "
              f"through degree {degree}", flush=True)

    basis, processed, terminated = complete(
        relations, key, args.max_degree, progress=progress)

    leads = sorted(basis, key=key)
    by_degree = {}
    for lead in leads:
        by_degree.setdefault(len(lead), 0)
        by_degree[len(lead)] += 1
    print(f"order = deg-lex {args.order};  processed degree = {processed};"
          f"  terminated = {terminated}")
    print("basis leading words by degree: "
          + ", ".join(f"{d}:{c}" for d, c in sorted(by_degree.items())))

    if terminated:
        # rigorous full comparison
        probe = 3 + 2 * (len(leads) * 4 + 3)  # safe overshoot of 2(s+3)
        series, states = normal_word_counts(leads, probe)
        probe = 2 * (states + 3) + 3
        series, states = normal_word_counts(leads, probe)
        target = target_series(probe)
        agree = series == target
        first_bad = next(
            (n for n in range(probe + 1) if series[n] != target[n]), None)
        verdict = "STRONGLY_FREE" if agree else "NOT_STRONGLY_FREE"
        print(f"automaton states = {states}; comparison through degree "
              f"{probe} (rigorous bound 2(s+3))")
        if first_bad is not None:
            print(f"first deviation at degree {first_bad}: "
                  f"a_n = {series[first_bad]}, target {target[first_bad]}")
    else:
        # exact prefix: leading words of the partial basis determine
        # dim (A/I)_n for n <= processed degree
        series, states = normal_word_counts(leads, processed)
        target = target_series(processed)
        agree = series == target
        verdict = ("INCONCLUSIVE_SERIES_MATCHES" if agree
                   else "NOT_STRONGLY_FREE")
        print(f"series verified through degree {processed}: "
              f"{'matches target' if agree else 'DEVIATES'}")
    print(f"VERDICT: {verdict}")
    if args.json:
        print(json.dumps({
            "order": args.order,
            "processed_degree": processed,
            "terminated": terminated,
            "leading_words": leads,
            "verdict": verdict,
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
