#!/usr/bin/env python3
"""Head-word engine behind verify_flag_certificate.py (plain Python).

A target is a finite, factor-free list of words over the letters 1,2,3
(the intended head words).  The static check verifies that the number
of words of each degree with no factor in the list satisfies the
recurrence of 1/(1-3z+3z^3) (Lemma A.4 of the paper) and lists all
overlaps up to degree 2*maxlen-1.  The dynamic test reduces the three
cubic relations to distinct head words, runs the overlap completion
(Lemma A.5 / Algorithm A.6 of the paper) and checks that every overlap
above the largest head degree reduces to zero.

This module is a library: the completion engine and the finite-field
arithmetic used by tools/verify_flag_certificate.py.  It has no
command-line interface in the published package.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ALPHABET = "123"
W0 = (
    "111", "112", "113", "1213", "1223", "1232", "12123",
    "12213", "12223", "122123", "1222123", "12221323",
)


def word_key(word: str):
    return len(word), word


class FiniteField:
    """F_p[t]/(m) with the coefficient-list convention of the certificates.

    Modulus coefficients are listed from the constant term upwards; the
    polynomials are the ones documented in FORMAT.md, keyed by (p, degree).
    """

    MODULI = {
        (3, 1): (0, 1),
        (3, 2): (2, 2, 1),             # t^2+2t+2
        (3, 3): (1, 2, 0, 1),          # t^3+2t+1
        (3, 4): (2, 1, 0, 0, 1),       # t^4+t+2
        (3, 5): (1, 2, 0, 0, 0, 1),        # t^5+2t+1
        (3, 6): (2, 1, 0, 0, 0, 0, 1),     # t^6+t+2
        (3, 7): (1, 0, 2, 0, 0, 0, 0, 1),  # t^7+2t^2+1
        (3, 8): (2, 0, 0, 1, 0, 0, 0, 0, 1),  # t^8+t^3+2
        (5, 1): (0, 1),
        (5, 2): (3, 0, 1),             # t^2+3 (2 is not a square mod 5)
    }

    def __init__(self, degree: int, prime: int = 3):
        if (prime, degree) not in self.MODULI:
            raise ValueError(f"no documented modulus for p={prime}, degree {degree}")
        self.p = prime
        self.degree = degree
        self.modulus = self.MODULI[(prime, degree)]
        self.size = prime ** degree
        self.zero = 0
        self.one = 1
        self.elements = tuple(range(self.size))
        self._digits = [
            tuple((value // (prime ** i)) % prime for i in range(degree))
            for value in self.elements
        ]
        self._neg = [
            self._encode(tuple(-a % prime for a in self._digits[x]))
            for x in self.elements
        ]
        self._pow = [prime ** i for i in range(degree)]
        if self.size <= 81:
            # small fields: full tables
            self._add = [
                [self._encode(tuple((a + b) % prime for a, b in
                                    zip(self._digits[x], self._digits[y])))
                 for y in self.elements]
                for x in self.elements
            ]
            self._mul = [
                [self._multiply(x, y) for y in self.elements]
                for x in self.elements
            ]
            self._inv = [None] * self.size
            for x in self.elements[1:]:
                self._inv[x] = next(
                    (y for y in self.elements if self._mul[x][y] == self.one),
                    None,
                )
                if self._inv[x] is None:
                    raise ValueError(
                        f"modulus for p={prime}, degree {degree} is not irreducible"
                    )
            self._log = None
        else:
            # larger fields: logarithm tables for a primitive element;
            # addition digit by digit.
            order = self.size - 1
            self._add = None
            self._mul = None
            self._inv = None
            for g in self.elements[2:]:
                exp = [self.one]
                x = g
                while x != self.one and len(exp) < order:
                    exp.append(x)
                    x = self._multiply(x, g)
                if len(exp) == order and x == self.one:
                    break
            else:
                raise ValueError(f"no primitive element for p={prime}, degree {degree}")
            self._exp = exp + exp          # doubled, saves the reduction mod order
            self._log = [None] * self.size
            for i, value in enumerate(exp):
                self._log[value] = i
            if any(self._log[v] is None for v in self.elements[1:]):
                raise ValueError(f"modulus for p={prime}, degree {degree} is not irreducible")
            self._order = order

    def _encode(self, coefficients):
        return sum(c * self.p ** i for i, c in enumerate(coefficients))

    def _multiply(self, x, y):
        degree, prime = self.degree, self.p
        product = [0] * (2 * degree - 1)
        for i, a in enumerate(self._digits[x]):
            for j, b in enumerate(self._digits[y]):
                product[i + j] = (product[i + j] + a * b) % prime
        for k in range(len(product) - 1, degree - 1, -1):
            coefficient = product[k]
            if coefficient:
                for i in range(degree):
                    product[k - degree + i] = (
                        product[k - degree + i]
                        - coefficient * self.modulus[i]
                    ) % prime
        return self._encode(tuple(product[:degree]))

    def elt(self, value):
        if isinstance(value, int):
            return value % self.p
        if not isinstance(value, (list, tuple)) or len(value) != self.degree:
            raise ValueError(f"not an element of F_({self.p}^{self.degree}): {value!r}")
        return self._encode(tuple(int(c) % self.p for c in value))

    def add(self, x, y):
        if self._add is not None:
            return self._add[x][y]
        dx, dy = self._digits[x], self._digits[y]
        return sum(((a + b) % self.p) * w for a, b, w in zip(dx, dy, self._pow))

    def neg(self, x):
        return self._neg[x]

    def sub(self, x, y):
        return self.add(x, self._neg[y])

    def mul(self, x, y):
        if self._mul is not None:
            return self._mul[x][y]
        if x == 0 or y == 0:
            return 0
        return self._exp[self._log[x] + self._log[y]]

    def inv(self, x):
        if x == self.zero:
            raise ZeroDivisionError
        if self._inv is not None:
            return self._inv[x]
        return self._exp[(self._order - self._log[x]) % self._order]

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    def is_zero(self, x):
        return x == self.zero

    def show(self, x):
        return x if self.degree == 1 else list(self._digits[x])

    def frobenius(self, x):
        result = self.one
        for _ in range(self.p):
            result = self.mul(result, x)
        return result


def infer_degree(raw_basis) -> int:
    sample = raw_basis[0][0]
    return len(sample) if isinstance(sample, (list, tuple)) else 1


def clean(poly, field: FiniteField):
    return {word: value for word, value in poly.items()
            if not field.is_zero(value)}


def sub_scaled(target, factor, source, field: FiniteField):
    if field.is_zero(factor):
        return
    for word, coefficient in source.items():
        value = field.sub(target.get(word, field.zero),
                          field.mul(factor, coefficient))
        if field.is_zero(value):
            target.pop(word, None)
        else:
            target[word] = value


def first_divisor(word, basis):
    for lead in basis:
        position = word.find(lead)
        if position >= 0:
            return position, lead
    return None


def reduce_poly(poly, basis, field: FiniteField):
    """Homogene Links-rechts-Reduktion in Dp (Grad, dann 1>2>3)."""
    poly = dict(poly)
    queue = list(poly)
    heapq.heapify(queue)
    normal = set()
    while queue:
        word = heapq.heappop(queue)
        if word not in poly or word in normal:
            continue
        hit = first_divisor(word, basis)
        if hit is None:
            normal.add(word)
            continue
        position, lead = hit
        left = word[:position]
        right = word[position + len(lead):]
        factor = poly[word]
        shifted = {left + w + right: c for w, c in basis[lead].items()}
        sub_scaled(poly, factor, shifted, field)
        for fresh in shifted:
            if fresh in poly and fresh not in normal:
                heapq.heappush(queue, fresh)
    return poly


def monic_at(poly, lead, field: FiniteField):
    inverse = field.inv(poly[lead])
    return {word: field.mul(value, inverse)
            for word, value in poly.items()}


def determinant3(matrix, field: FiniteField):
    a = matrix
    positive = field.add(
        field.add(field.mul(a[0][0], field.mul(a[1][1], a[2][2])),
                  field.mul(a[0][1], field.mul(a[1][2], a[2][0]))),
        field.mul(a[0][2], field.mul(a[1][0], a[2][1])))
    negative = field.add(
        field.add(field.mul(a[0][2], field.mul(a[1][1], a[2][0])),
                  field.mul(a[0][1], field.mul(a[1][0], a[2][2]))),
        field.mul(a[0][0], field.mul(a[1][2], a[2][1])))
    return field.sub(positive, negative)


def transform_tensor(tensor, raw_basis, field: FiniteField):
    g = [[field.elt(x) for x in row] for row in raw_basis]
    transformed = []
    for relation in range(3):
        row = []
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    total = field.zero
                    for i in range(3):
                        for j in range(3):
                            for k in range(3):
                                factor = field.mul(
                                    g[i][a], field.mul(g[j][b], g[k][c])
                                )
                                total = field.add(
                                    total,
                                    field.mul(tensor[relation][9*i+3*j+k],
                                              factor),
                                )
                    row.append(total)
        transformed.append(row)
    return transformed


def rows_to_polynomials(matrix, field: FiniteField):
    words = ["".join(word) for word in itertools.product(ALPHABET, repeat=3)]
    return [clean(dict(zip(words, row)), field) for row in matrix]


def canonical_relations(relations, field: FiniteField):
    """Volle RREF; ihre Pivotwoerter sind die eindeutigen Grad-3-Koepfe."""
    rows = [dict(row) for row in relations]
    words = ["".join(word) for word in itertools.product(ALPHABET, repeat=3)]
    pivots = []
    pivot_row = 0
    for word in words:
        found = next((i for i in range(pivot_row, len(rows))
                      if not field.is_zero(rows[i].get(word, field.zero))), None)
        if found is None:
            continue
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
        rows[pivot_row] = monic_at(rows[pivot_row], word, field)
        for i in range(len(rows)):
            if i != pivot_row:
                sub_scaled(rows[i], rows[i].get(word, field.zero),
                           rows[pivot_row], field)
        pivots.append(word)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    basis = {pivot: rows[i] for i, pivot in enumerate(pivots)}
    return tuple(pivots), basis


def pivot_minor(relations, pivots, field: FiniteField):
    matrix = [[relations[i].get(word, field.zero) for word in pivots]
              for i in range(3)]
    return determinant3(matrix, field)


def composition(left, right, composite, basis, field: FiniteField):
    suffix = composite[len(left):]
    prefix = composite[:len(composite) - len(right)]
    result = {word + suffix: coefficient
              for word, coefficient in basis[left].items()}
    shifted = {prefix + word: coefficient
               for word, coefficient in basis[right].items()}
    sub_scaled(result, field.one, shifted, field)
    return reduce_poly(result, basis, field)


def overlap_composites(left, right):
    for length in range(1, min(len(left), len(right))):
        if left[-length:] == right[:length]:
            yield length, left + right[length:]


def ambiguity_rows(basis, degree, field: FiniteField):
    rows = []
    leads = sorted(basis, key=word_key)
    for left in leads:
        for right in leads:
            for _, composite in overlap_composites(left, right):
                if len(composite) == degree:
                    rows.append(composition(left, right, composite,
                                            basis, field))
    return rows


def insert_row_space(rows, basis, field: FiniteField):
    new_heads = []
    for row in rows:
        row = reduce_poly(row, basis, field)
        if not row:
            continue
        lead = min(row)
        basis[lead] = monic_at(row, lead, field)
        new_heads.append(lead)
    return tuple(sorted(new_heads))


def canonical_words(raw_words):
    if not isinstance(raw_words, (list, tuple)):
        raise ValueError("W muss eine JSON-Liste sein")
    if len(set(raw_words)) != len(raw_words):
        raise ValueError("W enthaelt Duplikate")
    words = tuple(sorted(raw_words, key=word_key))
    if not words:
        raise ValueError("W darf nicht leer sein")
    if any(not isinstance(word, str) or not word
           or any(letter not in ALPHABET for letter in word)
           for word in words):
        raise ValueError("Kopfwoerter duerfen nur die Ziffern 1,2,3 enthalten")
    cubic = tuple(word for word in words if len(word) == 3)
    if len(cubic) != 3 or any(len(word) < 3 for word in words):
        raise ValueError("a target needs exactly three head words of degree 3")
    for i, short in enumerate(words):
        for j, long in enumerate(words):
            if i != j and short in long:
                raise ValueError(f"W ist nicht faktor-minimal: {short} | {long}")
    return words


def target_coefficients(limit):
    values = [1, 3, 9]
    while len(values) <= limit:
        values.append(3 * values[-1] - 3 * values[-3])
    return values[:limit + 1]


def avoidance_proof(words):
    prefixes = {""}
    for word in words:
        prefixes.update(word[:i] for i in range(1, len(word)))
    states = tuple(sorted(prefixes, key=word_key))
    state_index = {state: i for i, state in enumerate(states)}
    forbidden = set(words)

    def transition(state, letter):
        extended = state + letter
        if any(extended.endswith(word) for word in forbidden):
            return None
        while extended not in prefixes:
            extended = extended[1:]
        return state_index[extended]

    transitions = [[transition(state, letter) for letter in ALPHABET]
                   for state in states]
    # Wie beim W0-Nachweis wird inklusive n=s+3 geprueft.
    window_end = len(states) + 3
    counts = [1]
    distribution = [0] * len(states)
    distribution[state_index[""]] = 1
    for _ in range(window_end):
        fresh = [0] * len(states)
        for source, count in enumerate(distribution):
            if count:
                for target in transitions[source]:
                    if target is not None:
                        fresh[target] += count
        distribution = fresh
        counts.append(sum(distribution))
    target = target_coefficients(window_end)
    if counts != target:
        mismatch = next(i for i, pair in enumerate(zip(counts, target))
                        if pair[0] != pair[1])
        raise ValueError(
            f"falsche Vermeidungsreihe: erste Abweichung in Grad {mismatch} "
            f"({counts[mismatch]} != {target[mismatch]})"
        )
    return {
        "exact": True,
        "method": "prefix-automaton-plus-Cayley-Hamilton",
        "automaton_states": len(states),
        "states": list(states),
        "transitions_123": transitions,
        "recurrence": "a_n=3*a_(n-1)-3*a_(n-3)",
        "initial_values": [1, 3, 9],
        "window": [0, window_end],
        "coefficients": counts,
    }


def critical_composition_proof(words):
    rows = []
    for left in words:
        for right in words:
            for overlap, composite in overlap_composites(left, right):
                rows.append({
                    "left": left,
                    "right": right,
                    "overlap": overlap,
                    "composite": composite,
                    "degree": len(composite),
                })
    rows.sort(key=lambda row: (row["degree"], row["composite"],
                               row["left"], row["right"], row["overlap"]))
    maximum = max(map(len, words))
    bound = 2 * maximum - 1
    if any(row["degree"] > bound for row in rows):
        raise AssertionError("kritisches Kompositum oberhalb der Schranke")
    return {
        "exact": True,
        "method": "all-proper-suffix-prefix-overlaps",
        "highest_head_degree": maximum,
        "diamond_bound": bound,
        "count": len(rows),
        "maximum_composite_degree": max(
            (row["degree"] for row in rows), default=0
        ),
        "above_highest_head_degree": [
            row for row in rows if row["degree"] > maximum
        ],
        "all": rows,
    }


def derived_profiles(words):
    maximum = max(map(len, words))
    profiles = {}
    earlier = []
    for degree in range(3, maximum + 1):
        expected = tuple(word for word in words if len(word) == degree)
        normal = [
            "".join(word)
            for word in itertools.product(ALPHABET, repeat=degree)
            if not any(head in "".join(word) for head in earlier)
        ]
        positions = [normal.index(word) + 1 for word in expected]
        ambiguity_count = 0
        for left in earlier:
            for right in earlier:
                ambiguity_count += sum(
                    len(composite) == degree
                    for _, composite in overlap_composites(left, right)
                )
        profiles[str(degree)] = {
            "heads": list(expected),
            "normal_columns_before_insertion": len(normal),
            "pivot_positions": positions,
            "critical_rows": ambiguity_count if degree > 3 else 0,
        }
        earlier.extend(expected)
    return profiles


def precheck_cell(raw_words):
    words = canonical_words(raw_words)
    return {
        "head_words": list(words),
        "avoidance_series": avoidance_proof(words),
        "pivot_profiles": derived_profiles(words),
        "critical_compositions": critical_composition_proof(words),
    }


def prepare_cell(raw, default_id=None):
    if isinstance(raw, dict):
        words = raw.get("head_words", raw.get("kopfwoerter", raw.get("words")))
        cell_id = raw.get("id", default_id)
    else:
        words = raw
        cell_id = default_id
    words = canonical_words(words)
    proof = precheck_cell(words)
    expected = {
        degree: tuple(word for word in words if len(word) == degree)
        for degree in range(3, max(map(len, words)) + 1)
    }
    return {
        "id": cell_id,
        "head_words": words,
        "initial_heads": expected[3],
        "expected": expected,
        "max_degree": max(map(len, words)),
        "proof": proof,
    }


def _failure(cell, stage, **extra):
    result = {
        "W_regular": False,
        "in_cell": False,
        "cell_id": cell.get("id"),
        "failed_stage": stage,
        "failed": stage,
    }
    if cell["head_words"] == W0:
        result["in_W0"] = False
    result.update(extra)
    return result


def test_transformed(transformed, field: FiniteField, raw_cell,
                     details=False, check_termination=True):
    cell = raw_cell if isinstance(raw_cell, dict) and "expected" in raw_cell \
        else prepare_cell(raw_cell)
    relations = rows_to_polynomials(transformed, field)
    delta = pivot_minor(relations, cell["initial_heads"], field)
    values = {"P3_minor": field.show(delta)}
    if cell["head_words"] == W0:
        values["P0_Delta"] = field.show(delta)
    if field.is_zero(delta):
        stage = "P0_Delta" if cell["head_words"] == W0 else "P3_minor"
        return _failure(cell, stage, values=values)

    pivots, basis = canonical_relations(relations, field)
    if pivots != cell["initial_heads"]:
        return _failure(cell, "R3", expected=list(cell["initial_heads"]),
                        found=list(pivots), values=values)

    rank_profile = {
        "3": {"rank": len(pivots), "pivots": list(pivots)}
    }
    for degree in range(4, cell["max_degree"] + 1):
        rows = ambiguity_rows(basis, degree, field)
        found = insert_row_space(rows, basis, field)
        expected = cell["expected"].get(degree, ())
        rank_profile[str(degree)] = {
            "ambiguity_rows": len(rows),
            "rank": len(found),
            "pivots": list(found),
        }
        if found != expected:
            return _failure(cell, f"R{degree}", expected=list(expected),
                            found=list(found), values=values,
                            rank_profile=rank_profile)

    termination = []
    if check_termination:
        bound = 2 * cell["max_degree"] - 1
        for degree in range(cell["max_degree"] + 1, bound + 1):
            rows = ambiguity_rows(basis, degree, field)
            found = insert_row_space(rows, basis, field)
            termination.append({
                "degree": degree,
                "critical_rows": len(rows),
                "nonzero_remainders": len(found),
                "pivots": list(found),
            })
            if found:
                return _failure(cell, f"T{degree}", expected=[],
                                found=list(found), values=values,
                                rank_profile=rank_profile,
                                termination=termination)

    result = {
        "W_regular": True,
        "in_cell": True,
        "cell_id": cell.get("id"),
        "failed_stage": None,
        "failed": None,
        "highest_head_degree": cell["max_degree"],
    }
    if cell["head_words"] == W0:
        result["in_W0"] = True
    if details:
        result.update(
            values=values,
            head_words=list(cell["head_words"]),
            rank_profile=rank_profile,
            termination=termination,
            precheck=cell["proof"],
        )
    return result


def test_cells_transformed(transformed, field: FiniteField, raw_cells,
                           details=False):
    """Gemeinsamer Praefixtrace; ein Grad wird je ueberlebendem Ast berechnet."""
    cells = [cell if isinstance(cell, dict) and "expected" in cell
             else prepare_cell(cell, f"cell-{index + 1}")
             for index, cell in enumerate(raw_cells)]
    relations = rows_to_polynomials(transformed, field)
    pivots, initial_basis = canonical_relations(relations, field)
    candidates = [cell for cell in cells if cell["initial_heads"] == pivots]
    trace = {"3": {"rank": len(pivots), "pivots": list(pivots),
                   "candidates": [cell["id"] for cell in candidates]}}
    if not candidates:
        return {"found": False, "cell_id": None, "failed_stage": "R3",
                "actual_pivots": list(pivots), "trace": trace if details else None}

    # Kandidaten mit demselben Grad-3-Praefix teilen exakt dieselbe Basis.
    basis = dict(initial_basis)
    maximum = max(cell["max_degree"] for cell in candidates)
    for degree in range(4, maximum + 1):
        rows = ambiguity_rows(basis, degree, field)
        found = insert_row_space(rows, basis, field)
        candidates = [
            cell for cell in candidates
            if cell["expected"].get(degree, ()) == found
        ]
        trace[str(degree)] = {
            "ambiguity_rows": len(rows), "pivots": list(found),
            "candidates": [cell["id"] for cell in candidates],
        }
        if not candidates:
            return {"found": False, "cell_id": None,
                    "failed_stage": f"R{degree}",
                    "trace": trace if details else None}

    # Verschiedene IDs duerfen dasselbe W benennen; Prioritaet ist Listenfolge.
    cell = candidates[0]
    bound = 2 * cell["max_degree"] - 1
    termination = []
    for degree in range(maximum + 1, bound + 1):
        rows = ambiguity_rows(basis, degree, field)
        found = insert_row_space(rows, basis, field)
        termination.append({"degree": degree, "critical_rows": len(rows),
                            "pivots": list(found)})
        if found:
            return {"found": False, "cell_id": None,
                    "failed_stage": f"T{degree}",
                    "trace": trace if details else None,
                    "termination": termination if details else None}
    result = {
        "found": True,
        "cell_id": cell["id"],
        "head_words": list(cell["head_words"]),
        "highest_head_degree": cell["max_degree"],
        "failed_stage": None,
    }
    if details:
        result.update(trace=trace, termination=termination)
    return result


def load_json_argument(argument):
    path_text = argument[1:] if argument.startswith("@") else argument
    path = Path(path_text)
    if path.exists():
        return json.loads(path.read_text()), path
    return json.loads(argument), None


def select_cell(data, cell_id=None):
    if isinstance(data, list):
        return data
    if not isinstance(data, dict):
        raise ValueError("invalid target JSON")
    if "cells" in data:
        if cell_id is None:
            raise ValueError("Katalogdatei verlangt --cell ID")
        selected = next((cell for cell in data["cells"]
                         if cell.get("id") == cell_id), None)
        if selected is None:
            raise KeyError(cell_id)
        return selected
    return data


def load_basis(argument):
    data, _ = load_json_argument(argument)
    return data["basis"] if isinstance(data, dict) and "basis" in data else data


def load_tensor(discriminant, extra_paths=()):
    paths = [Path(path) for path in extra_paths]
    paths.extend((
        HERE / "fahnensuche" / "tensoren-rest.json",
        HERE / "verifikation-block0" / "data" / "tensoren-block0.json",
    ))
    key = str(discriminant)
    for path in paths:
        if not path.exists():
            continue
        data = json.loads(path.read_text())
        if key in data:
            return data[key], path
    raise KeyError(f"Tensor {discriminant} in keiner Tensor-Datei gefunden")


# The command-line interface of the private tool is removed in this
# published copy: the module is a library used by verify_flag_certificate.py,
# which performs the schema and provenance checks.

if __name__ == "__main__":
    import sys
    print("flag_cell_checker.py is a library; use tools/verify_flag_certificate.py", file=sys.stderr)
    sys.exit(2)
