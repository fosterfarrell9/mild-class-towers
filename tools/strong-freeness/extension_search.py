#!/usr/bin/env python3
"""Exact Anick searches over F_3 and F_9 modulo the diagonal torus.

The matrix convention here is X_i -> sum_a M[i,a] Y_a.  Consequently
right multiplication by a diagonal matrix rescales output variables, and
projectively normalized columns give representatives of GL_3(F_q)/T.
"""

from __future__ import annotations

import json
import math
import random
import sys
import time
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)]

from admissible import (  # noqa: E402
    combinatorially_free,
    decode_word,
    gl3_matrices,
    leading_words,
    random_tensor,
    transform_tensor,
    word_name,
)

RESULTS = HERE / "results"
F3_SAMPLE_SEED = 20260804
F27_SAMPLE_SEED = 20260805


class FiniteField:
    """Small prime-power field F_3[t]/(modulus), with packed elements."""

    def __init__(self, modulus: tuple[int, ...]):
        if modulus[-1] % 3 != 1:
            raise ValueError("modulus must be monic")
        self.modulus = tuple(x % 3 for x in modulus)
        self.degree = len(modulus) - 1
        self.q = 3 ** self.degree
        self.coefficients = [self._unpack_raw(x) for x in range(self.q)]
        self.add = [[self._pack_raw(tuple((a + b) % 3 for a, b in zip(
            self.coefficients[x], self.coefficients[y])))
                     for y in range(self.q)] for x in range(self.q)]
        self.neg = [self._pack_raw(tuple(-a % 3 for a in self.coefficients[x]))
                    for x in range(self.q)]
        self.mul = [[self._multiply_raw(x, y) for y in range(self.q)]
                    for x in range(self.q)]
        self.inv = [0] * self.q
        for x in range(1, self.q):
            self.inv[x] = next(y for y in range(1, self.q)
                               if self.mul[x][y] == 1)

    def _unpack_raw(self, value: int) -> tuple[int, ...]:
        out = []
        for _ in range(self.degree):
            out.append(value % 3)
            value //= 3
        return tuple(out)

    @staticmethod
    def _pack_raw(coefficients: tuple[int, ...]) -> int:
        value = 0
        place = 1
        for coefficient in coefficients:
            value += (coefficient % 3) * place
            place *= 3
        return value

    def _multiply_raw(self, left: int, right: int) -> int:
        a = self.coefficients[left]
        b = self.coefficients[right]
        work = [0] * (2 * self.degree - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                work[i + j] = (work[i + j] + ai * bj) % 3
        for power in range(len(work) - 1, self.degree - 1, -1):
            coefficient = work[power] % 3
            if coefficient:
                for i in range(self.degree):
                    work[power - self.degree + i] = (
                        work[power - self.degree + i]
                        - coefficient * self.modulus[i]
                    ) % 3
        return self._pack_raw(tuple(work[:self.degree]))

    def element_name(self, value: int) -> str:
        coefficients = self.coefficients[value]
        terms = []
        for power, coefficient in enumerate(coefficients):
            if not coefficient:
                continue
            atom = "1" if power == 0 else "t" if power == 1 else f"t^{power}"
            terms.append(atom if coefficient == 1 else f"2*{atom}")
        return "+".join(terms) if terms else "0"


def determinant_columns(a, b, c, field: FiniteField) -> int:
    """Determinant of the matrix with columns a,b,c."""
    add, mul, neg = field.add, field.mul, field.neg
    positive = add[mul[a[0]][mul[b[1]][c[2]]]][
        add[mul[b[0]][mul[c[1]][a[2]]]][mul[c[0]][mul[a[1]][b[2]]]]]
    negative = add[mul[c[0]][mul[b[1]][a[2]]]][
        add[mul[b[0]][mul[a[1]][c[2]]]][mul[a[0]][mul[c[1]][b[2]]]]]
    return add[positive][neg[negative]]


def projective_vectors(field: FiniteField) -> list[tuple[int, int, int]]:
    vectors = []
    for vector in product(range(field.q), repeat=3):
        if not any(vector):
            continue
        first = next(x for x in vector if x)
        normalized = tuple(field.mul[field.inv[first]][x] for x in vector)
        if normalized == vector:
            vectors.append(vector)
    assert len(vectors) == field.q * field.q + field.q + 1
    return vectors


def coefficient_tables(tensor, vectors, field: FiniteField):
    """Tabulate each trilinear relation on all projective vector triples."""
    length = len(vectors)
    add, mul = field.add, field.mul
    tables = []
    for row in tensor:
        pair_values = [[[0, 0, 0] for _ in range(length)]
                       for _ in range(length)]
        for u_index, u in enumerate(vectors):
            for v_index, v in enumerate(vectors):
                output = pair_values[u_index][v_index]
                for k in range(3):
                    total = 0
                    for i in range(3):
                        for j in range(3):
                            coefficient = row[(i * 3 + j) * 3 + k]
                            if coefficient:
                                term = mul[mul[u[i]][v[j]]][coefficient]
                                total = add[total][term]
                    output[k] = total
        table = bytearray(length ** 3)
        offset = 0
        for u_index in range(length):
            for v_index in range(length):
                contracted = pair_values[u_index][v_index]
                for w in vectors:
                    total = add[mul[contracted[0]][w[0]]][
                        add[mul[contracted[1]][w[1]]][mul[contracted[2]][w[2]]]]
                    table[offset] = total
                    offset += 1
        tables.append(table)
    return tables


def projective_ids(field: FiniteField) -> list[int]:
    out = [0] * (field.q ** 3)
    for a, b, c in product(range(field.q), repeat=3):
        packed = (a * field.q + b) * field.q + c
        if not (a or b or c):
            out[packed] = -1
            continue
        first = a or b or c
        inverse = field.inv[first]
        normalized = tuple(field.mul[inverse][x] for x in (a, b, c))
        out[packed] = (normalized[0] * field.q + normalized[1]) * field.q + normalized[2]
    return out


def coset_anick_search(tensor, field: FiniteField, limit: int | None = None,
                       seed: int | None = None):
    """Search normalized-column representatives of GL_3(F_q)/T exactly.

    With ``limit`` this samples representatives uniformly with replacement;
    the exhaustive path is used for q=3 and q=9.
    """
    started = time.monotonic()
    vectors = projective_vectors(field)
    tables = coefficient_tables(tensor, vectors, field)
    table_seconds = time.monotonic() - started
    length = len(vectors)
    normalized_coefficients = projective_ids(field)
    descending_words = [(word, decode_word(word)) for word in range(26, -1, -1)]

    valid_thirds = {}
    for first in range(length):
        for second in range(length):
            if first == second:
                continue
            choices = [third for third in range(length)
                       if determinant_columns(vectors[first], vectors[second],
                                              vectors[third], field)]
            valid_thirds[first, second] = choices
    expected = (field.q * field.q + field.q + 1) * (
        field.q * field.q + field.q) * field.q * field.q
    total_representatives = sum(map(len, valid_thirds.values()))
    assert total_representatives == expected

    def check(first: int, second: int, third: int):
        columns = (first, second, third)
        selected_vectors = []
        leaders = []
        for word, (a, b, c) in descending_words:
            index = ((columns[a] * length + columns[b]) * length + columns[c])
            vector = (tables[0][index], tables[1][index], tables[2][index])
            if not any(vector):
                continue
            packed = (vector[0] * field.q + vector[1]) * field.q + vector[2]
            if not selected_vectors:
                selected_vectors.append(vector)
                leaders.append(word)
            elif len(selected_vectors) == 1:
                previous = selected_vectors[0]
                previous_packed = ((previous[0] * field.q + previous[1])
                                   * field.q + previous[2])
                if normalized_coefficients[packed] != normalized_coefficients[previous_packed]:
                    selected_vectors.append(vector)
                    leaders.append(word)
            elif determinant_columns(selected_vectors[0], selected_vectors[1],
                                     vector, field):
                leaders.append(word)
                assert len(leaders) == 3
                return leaders if combinatorially_free(leaders) else None
        raise AssertionError("invertible change lost relation rank")

    tested = 0
    witness = None
    witness_count = 0
    if limit is None:
        for first in range(length):
            for second in range(length):
                if first == second:
                    continue
                for third in valid_thirds[first, second]:
                    tested += 1
                    leaders = check(first, second, third)
                    if leaders is not None:
                        witness_count += 1
                        if witness is None:
                            witness = (first, second, third, leaders)
        exhaustive = True
    else:
        rng = random.Random(seed)
        pairs = list(valid_thirds)
        for _ in range(limit):
            first, second = pairs[rng.randrange(len(pairs))]
            choices = valid_thirds[first, second]
            third = choices[rng.randrange(len(choices))]
            tested += 1
            leaders = check(first, second, third)
            if leaders is not None:
                witness_count += 1
                if witness is None:
                    witness = (first, second, third, leaders)
        exhaustive = False

    elapsed = time.monotonic() - started
    output = {
        "field_order": field.q,
        "modulus_coefficients_low_to_high": field.modulus,
        "torus_order": (field.q - 1) ** 3,
        "gl3_order": (field.q ** 3 - 1) * (field.q ** 3 - field.q)
                     * (field.q ** 3 - field.q ** 2),
        "coset_representatives": total_representatives,
        "tested": tested,
        "exhaustive": exhaustive and tested == total_representatives,
        "found": witness is not None,
        "witness_representatives": witness_count,
        "precomputation_seconds": table_seconds,
        "elapsed_seconds": elapsed,
        "seconds_per_test_amortized": elapsed / tested,
        "seed": seed,
    }
    if witness:
        first, second, third, leaders = witness
        columns = [vectors[index] for index in (first, second, third)]
        matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
        output.update({
            "matrix_packed": matrix,
            "matrix_polynomial": [[field.element_name(x) for x in row]
                                  for row in matrix],
            "leaders": [word_name(word) for word in leaders],
        })
    else:
        output.update({"matrix_packed": None, "matrix_polynomial": None,
                       "leaders": []})
    return output


def sampled_coset_anick_search(tensor, field: FiniteField, samples: int,
                                seed: int):
    """Uniformly sample projective bases without a cubic-size table.

    For each ordered pair of distinct projective points, exactly q^2 points
    lie outside its projective line.  Rejection sampling the third column is
    therefore uniform on GL_3(F_q)/T.  Admissibility gives
    f(a,b,c)=f(c,b,a), so two quadratic lookup tables cover every word with
    a repeated column; only the six all-distinct patterns are evaluated
    directly.
    """
    started = time.monotonic()
    vectors = projective_vectors(field)
    length = len(vectors)
    add, mul = field.add, field.mul

    adjacent_tables = []  # A[a,b] = f(a,a,b)
    palindrome_tables = []  # B[a,b] = f(a,b,a)
    sparse_rows = []
    for row in tensor:
        sparse_rows.append([(i, j, k, coefficient)
                            for i, j, k in product(range(3), repeat=3)
                            if (coefficient := row[(i * 3 + j) * 3 + k])])
        adjacent = bytearray(length * length)
        palindrome = bytearray(length * length)
        for a_index, a in enumerate(vectors):
            contract_last = [0, 0, 0]
            contract_middle = [0, 0, 0]
            for i, j, k, coefficient in sparse_rows[-1]:
                aa_last = mul[mul[a[i]][a[j]]][coefficient]
                contract_last[k] = add[contract_last[k]][aa_last]
                aa_middle = mul[mul[a[i]][a[k]]][coefficient]
                contract_middle[j] = add[contract_middle[j]][aa_middle]
            offset = a_index * length
            for b_index, b in enumerate(vectors):
                adjacent[offset + b_index] = add[mul[contract_last[0]][b[0]]][
                    add[mul[contract_last[1]][b[1]]][mul[contract_last[2]][b[2]]]]
                palindrome[offset + b_index] = add[mul[contract_middle[0]][b[0]]][
                    add[mul[contract_middle[1]][b[1]]][mul[contract_middle[2]][b[2]]]]
        adjacent_tables.append(adjacent)
        palindrome_tables.append(palindrome)

    normalized_coefficients = projective_ids(field)
    descending_words = [(word, decode_word(word)) for word in range(26, -1, -1)]

    def direct_coefficient(relation: int, a_index: int, b_index: int,
                           c_index: int) -> int:
        a, b, c = vectors[a_index], vectors[b_index], vectors[c_index]
        total = 0
        for i, j, k, scalar in sparse_rows[relation]:
            total = add[total][mul[mul[mul[a[i]][b[j]]][c[k]]][scalar]]
        return total

    def coefficient(relation: int, a_index: int, b_index: int,
                    c_index: int) -> int:
        if a_index == b_index:
            return adjacent_tables[relation][a_index * length + c_index]
        if b_index == c_index:
            return adjacent_tables[relation][b_index * length + a_index]
        if a_index == c_index:
            return palindrome_tables[relation][a_index * length + b_index]
        return direct_coefficient(relation, a_index, b_index, c_index)

    # Deterministic table-vs-direct audit, covering all repetition patterns.
    for audit in range(24):
        a_index = (17 * audit + 1) % length
        b_index = (31 * audit + 2) % length
        c_index = (47 * audit + 3) % length
        triples = ((a_index, a_index, b_index),
                   (a_index, b_index, a_index),
                   (a_index, b_index, b_index),
                   (a_index, b_index, c_index))
        for indices in triples:
            for relation in range(3):
                if coefficient(relation, *indices) != direct_coefficient(
                        relation, *indices):
                    raise AssertionError("sample lookup-table audit failed")

    def check(columns):
        selected_vectors = []
        leaders = []
        for word, letters in descending_words:
            indices = tuple(columns[letter] for letter in letters)
            vector = tuple(coefficient(relation, *indices) for relation in range(3))
            if not any(vector):
                continue
            packed = (vector[0] * field.q + vector[1]) * field.q + vector[2]
            if not selected_vectors:
                selected_vectors.append(vector)
                leaders.append(word)
            elif len(selected_vectors) == 1:
                previous = selected_vectors[0]
                previous_packed = ((previous[0] * field.q + previous[1])
                                   * field.q + previous[2])
                if normalized_coefficients[packed] != normalized_coefficients[previous_packed]:
                    selected_vectors.append(vector)
                    leaders.append(word)
            elif determinant_columns(selected_vectors[0], selected_vectors[1],
                                     vector, field):
                leaders.append(word)
                return leaders if combinatorially_free(leaders) else None
        raise AssertionError("invertible change lost relation rank")

    precomputation_seconds = time.monotonic() - started
    rng = random.Random(seed)
    first_witness = None
    witness_count = 0
    rejected_third_draws = 0
    for tested in range(1, samples + 1):
        first = rng.randrange(length)
        second = rng.randrange(length - 1)
        if second >= first:
            second += 1
        while True:
            third = rng.randrange(length)
            if determinant_columns(vectors[first], vectors[second],
                                   vectors[third], field):
                break
            rejected_third_draws += 1
        leaders = check((first, second, third))
        if leaders is not None:
            witness_count += 1
            if first_witness is None:
                first_witness = (first, second, third, leaders, tested)
        if tested % 1_000_000 == 0:
            print(f"F27 sample {tested}/{samples}, witnesses={witness_count}",
                  flush=True)

    elapsed = time.monotonic() - started
    gl3_order = ((field.q ** 3 - 1) * (field.q ** 3 - field.q)
                 * (field.q ** 3 - field.q ** 2))
    cosets = gl3_order // ((field.q - 1) ** 3)
    output = {
        "field_order": field.q,
        "modulus_coefficients_low_to_high": field.modulus,
        "torus_order": (field.q - 1) ** 3,
        "gl3_order": gl3_order,
        "coset_representatives": cosets,
        "sampling": "uniform with replacement on projective bases",
        "tested": samples,
        "exhaustive": False,
        "seed": seed,
        "found": first_witness is not None,
        "witnesses_in_sample": witness_count,
        "observed_density": witness_count / samples,
        "zero_success_95_percent_upper_density": (
            1.0 - math.pow(0.05, 1.0 / samples) if witness_count == 0 else None
        ),
        "rejected_third_draws": rejected_third_draws,
        "precomputation_seconds": precomputation_seconds,
        "elapsed_seconds": elapsed,
        "seconds_per_test_amortized": elapsed / samples,
    }
    if first_witness:
        first, second, third, leaders, position = first_witness
        columns = [vectors[index] for index in (first, second, third)]
        matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
        output.update({
            "first_witness_sample": position,
            "matrix_packed": matrix,
            "matrix_polynomial": [[field.element_name(x) for x in row]
                                  for row in matrix],
            "leaders": [word_name(word) for word in leaders],
        })
    else:
        output.update({"first_witness_sample": None, "matrix_packed": None,
                       "matrix_polynomial": None, "leaders": []})
    return output


def full_f3_anick_search(tensor):
    """Enumerate all 11232 matrices, retaining the first and total witnesses."""
    started = time.monotonic()
    first = None
    witness_count = 0
    tested = 0
    for matrix in gl3_matrices():
        tested += 1
        leaders = leading_words(transform_tensor(tensor, matrix))
        if combinatorially_free(leaders):
            witness_count += 1
            if first is None:
                first = {"matrix": matrix,
                         "leaders": [word_name(word) for word in leaders],
                         "position": tested}
    assert tested == 11232
    return {"found": first is not None, "tested": tested, "exhaustive": True,
            "witness_matrices": witness_count, "first_witness": first,
            "elapsed_seconds": time.monotonic() - started}


def f3_sanity(pilot_tensor):
    rng = random.Random(F3_SAMPLE_SEED)
    samples = [random_tensor(rng) for _ in range(3)]
    cases = [("finite-random-0-known-strong", samples[0]),
             ("finite-random-1", samples[1]),
             ("finite-random-2", samples[2]),
             ("arithmetic-pilot", pilot_tensor)]
    field = FiniteField((0, 1))
    records = []
    for label, tensor in cases:
        full = full_f3_anick_search(tensor)
        coset = coset_anick_search(tensor, field)
        if full["found"] != coset["found"]:
            raise AssertionError(f"torus reduction disagrees for {label}")
        if label.endswith("known-strong") and not full["found"]:
            raise AssertionError("recorded known strong sample lost its witness")
        records.append({"case": label, "full_gl3": full, "torus_cosets": coset,
                        "agree_on_existence": True})
    return {
        "seed": F3_SAMPLE_SEED,
        "full_search_size": 11232,
        "coset_search_size": 1404,
        "cases": records,
        "all_agree": True,
    }


def main():
    RESULTS.mkdir(exist_ok=True)
    pilot_tensor = json.loads((RESULTS / "tensor.json").read_text())["tensor_3_by_27"]

    sanity = f3_sanity(pilot_tensor)
    (RESULTS / "f3-torus-sanity.json").write_text(
        json.dumps(sanity, indent=2) + "\n")

    f9 = FiniteField((1, 0, 1))  # t^2+1
    result = coset_anick_search(pilot_tensor, f9)
    result["interpretation"] = (
        "A witness proves strong freeness of the reconstructed tensor over F3 "
        "by scalar-extension invariance. It says nothing about the tower until "
        "the p=3 bridge theorem, including restricted cubes, is proved."
        if result["found"] else
        "No F9 Anick witness leaves strong freeness undecided; it is not a "
        "negative verdict."
    )
    (RESULTS / "f9-witness.json").write_text(json.dumps(result, indent=2) + "\n")

    f27_result = None
    if not result["found"]:
        f27 = FiniteField((2, 2, 0, 1))  # t^3-t-1
        f27_result = sampled_coset_anick_search(
            pilot_tensor, f27, samples=10_000_000, seed=F27_SAMPLE_SEED)
        f27_result["interpretation"] = (
            "A sampled witness is positive evidence and a certificate. No "
            "sampled witness is only a density bound, never a negative verdict."
        )
        (RESULTS / "f27-sample.json").write_text(
            json.dumps(f27_result, indent=2) + "\n")

    seeds = json.loads((RESULTS / "seeds.json").read_text())
    seeds.update({
        "f3_torus_sanity_seed": F3_SAMPLE_SEED,
        "f27_coset_sample_seed": F27_SAMPLE_SEED,
        "note": ("F3 and F9 searches are exhaustive and deterministic. "
                 "The F27 seed is used only if F9 has no witness."),
    })
    seeds["randomness_used"] = True
    seeds["seeds"] = ([F3_SAMPLE_SEED, F27_SAMPLE_SEED]
                      if f27_result is not None else [F3_SAMPLE_SEED])
    (RESULTS / "seeds.json").write_text(json.dumps(seeds, indent=2) + "\n")

    print(json.dumps({"f3_sanity": sanity["all_agree"], "f9": result,
                      "f27": f27_result}, indent=2))


if __name__ == "__main__":
    main()
