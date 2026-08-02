#!/usr/bin/env python3
"""Parameterized noncommutative Groebner test for strong freeness.

This is a fresh, local characteristic-independent port of the p=5 blueprint.
The coefficient prime is explicit (default 3); the p=5 production file is not
imported or modified.  For three cubic relations the strongly-free target is

    1 / (1 - 3 z + 3 z^3).

A Hilbert-series deviation is a rigorous negative certificate.  Terminating
completion followed by the finite automaton comparison is a positive or
negative decision.  A degree cutoff without deviation is inconclusive.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

LETTERS = "xyz"
PRIME = 3


def set_prime(prime: int) -> None:
    global PRIME
    if prime < 3 or any(prime % q == 0 for q in range(2, int(prime**0.5)+1)):
        raise ValueError("prime must be an odd prime")
    PRIME = prime


def poly_normalize(poly: dict) -> dict:
    return {w: c % PRIME for w, c in poly.items() if c % PRIME}


def poly_sub_scaled(target: dict, factor: int, source: dict) -> None:
    for word, coefficient in source.items():
        value = (target.get(word, 0) - factor * coefficient) % PRIME
        if value:
            target[word] = value
        else:
            target.pop(word, None)


def make_order_key(letter_order: str):
    rank = {letter: i for i, letter in enumerate(letter_order)}

    def key(word: str):
        return len(word), [rank[letter] for letter in word]
    return key


def leading_word(poly: dict, key) -> str:
    return max(poly, key=key)


def find_division(word: str, leads: dict):
    for lead in leads:
        position = word.find(lead)
        if position >= 0:
            return position, lead
    return None


def reduce_poly(poly: dict, basis: dict, key) -> dict:
    import heapq

    def heap_entry(word):
        degree, lex = key(word)
        return -degree, [-v for v in lex], word

    poly = dict(poly)
    heap = [heap_entry(word) for word in poly]
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
        left, right = word[:position], word[position+len(lead):]
        factor = poly[word]
        shifted = {left+w+right: c for w, c in basis[lead].items()}
        poly_sub_scaled(poly, factor, shifted)
        for fresh in shifted:
            if fresh in poly and fresh not in normal:
                heapq.heappush(heap, heap_entry(fresh))
    return poly


def monic(poly: dict, key) -> dict:
    lead = leading_word(poly, key)
    inverse = pow(poly[lead], PRIME - 2, PRIME)
    return {w: c*inverse % PRIME for w, c in poly.items()}


def overlap_ambiguities(left: str, right: str):
    for length in range(1, min(len(left), len(right))):
        if left[-length:] == right[:length]:
            yield left + right[length:]


def complete(relations, key, max_degree, progress=None):
    basis = {}
    pending = [poly_normalize(relation) for relation in relations]

    def insert_all():
        nonlocal pending
        changed = True
        while changed:
            changed = False
            queue, pending = pending, []
            for poly in queue:
                reduced = reduce_poly(poly, basis, key)
                if not reduced:
                    continue
                reduced = monic(reduced, key)
                lead = leading_word(reduced, key)
                retired = [old for old in basis if old != lead and lead in old]
                basis[lead] = reduced
                for old in retired:
                    pending.append(basis.pop(old))
                changed = True

    insert_all()
    degree = 3
    while degree < max_degree:
        degree += 1
        new_elements = []
        leads = sorted(basis, key=key)
        for left in leads:
            for right in leads:
                for composite in overlap_ambiguities(left, right):
                    if len(composite) != degree:
                        continue
                    g1, g2 = basis.get(left), basis.get(right)
                    if g1 is None or g2 is None:
                        continue
                    suffix = composite[len(left):]
                    prefix = composite[:len(composite)-len(right)]
                    s_poly = {w+suffix: c for w, c in g1.items()}
                    poly_sub_scaled(s_poly, 1, {prefix+w: c for w, c in g2.items()})
                    remainder = reduce_poly(s_poly, basis, key)
                    if remainder:
                        new_elements.append(remainder)
        if new_elements:
            pending.extend(new_elements)
            insert_all()
        if progress:
            progress(degree, basis)
        largest = max(len(word) for word in basis)
        if degree >= 2*largest - 1:
            return basis, degree, True
    return basis, max_degree, False


def build_automaton(leads):
    prefixes = {""}
    for lead in leads:
        prefixes.update(lead[:i] for i in range(1, len(lead)))
    states = sorted(prefix for prefix in prefixes
                    if find_division(prefix, dict.fromkeys(leads)) is None)
    index = {state: i for i, state in enumerate(states)}
    transitions = []
    for state in states:
        row = []
        for letter in LETTERS:
            word = state + letter
            if any(word.endswith(lead) for lead in leads):
                row.append(None)
            else:
                row.append(next(index[word[start:]]
                                for start in range(len(word)+1)
                                if word[start:] in index))
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
            for target in transitions[state]:
                if target is not None:
                    fresh[target] += value
        counts = fresh
        series.append(sum(counts))
    return series, len(states)


def target_series(length):
    values = [1, 3, 9]
    for degree in range(3, length+1):
        values.append(3*values[degree-1] - 3*values[degree-3])
    return values[:length+1]


def rows_to_polynomials(matrix):
    words = [LETTERS[i]+LETTERS[j]+LETTERS[k]
             for i in range(3) for j in range(3) for k in range(3)]
    return [{words[q]: value % PRIME for q, value in enumerate(row)
             if value % PRIME} for row in matrix]


def hilbert_verdict(matrix, max_degree=10, order="xyz"):
    key = make_order_key(order)
    basis, processed, terminated = complete(
        rows_to_polynomials(matrix), key, max_degree)
    leads = list(basis)
    series, states = normal_word_counts(leads, max(processed, 1))
    target = target_series(max(processed, 1))
    for degree in range(processed+1):
        if series[degree] != target[degree]:
            return {
                "verdict": "NOT_STRONGLY_FREE",
                "deviation_degree": degree,
                "processed_degree": processed,
                "terminated": terminated,
            }
    if terminated:
        bound = 2*(states+3)+3
        series, _ = normal_word_counts(leads, bound)
        target = target_series(bound)
        if series == target:
            return {"verdict": "STRONGLY_FREE", "deviation_degree": None,
                    "processed_degree": processed, "terminated": True}
        degree = next(i for i, (a, b) in enumerate(zip(series, target)) if a != b)
        return {"verdict": "NOT_STRONGLY_FREE", "deviation_degree": degree,
                "processed_degree": processed, "terminated": True}
    return {"verdict": "INCONCLUSIVE", "deviation_degree": None,
            "processed_degree": processed, "terminated": False}


def parse_matrix(text: str):
    rows = text.strip().lstrip("[").rstrip("]").split(";")
    matrix = [[int(v) % PRIME for v in row.split(",")] for row in rows]
    if len(matrix) != 3 or any(len(row) != 27 for row in matrix):
        raise ValueError("matrix must be 3 x 27")
    return matrix


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=3)
    parser.add_argument("--matrix")
    parser.add_argument("--matrix-file", type=Path)
    parser.add_argument("--max-degree", type=int, default=12)
    parser.add_argument("--order", default="xyz")
    args = parser.parse_args(argv)
    set_prime(args.prime)
    if bool(args.matrix) == bool(args.matrix_file):
        parser.error("give exactly one of --matrix and --matrix-file")
    matrix = parse_matrix(args.matrix if args.matrix else args.matrix_file.read_text())
    print(json.dumps(hilbert_verdict(matrix, args.max_degree, args.order)))


if __name__ == "__main__":
    main()
