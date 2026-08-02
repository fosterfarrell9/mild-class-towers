#!/usr/bin/env python3
"""Reconstruction of characteristic-three admissible cubic tensors.

This script determines minimal sets of projective D-evaluations by exact
linear algebra and verifies the formulas on seeded random tensors.  It writes
``results/reconstruction.json``.
"""

from __future__ import annotations

from itertools import combinations
import json
import random
from pathlib import Path

from admissible import (
    P, admissible_basis, contract_D, contract_row, diagonal_matrix,
    projective_points, random_tensor, rank, rref, validate_tensor, word_idx,
)

HERE = Path(__file__).resolve().parent
SEED = 20260802
BASIS_POINTS = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
PAIR_POINTS = [(1, 1, 0), (1, 0, 1), (0, 1, 1)]


def evaluation_matrix(points):
    """Matrix of row-tensor -> (D_x)_m for one relation coordinate."""
    basis = admissible_basis()
    return [[contract_row(row, point)[m] for row in basis]
            for point in points for m in range(3)]


def evaluation_rank(points) -> int:
    return rank(evaluation_matrix(points))


def minimal_families():
    points = projective_points()
    for size in range(1, len(points) + 1):
        full = [combo for combo in combinations(range(13), size)
                if evaluation_rank([points[i] for i in combo]) == 11]
        if full:
            return size, full
    raise AssertionError("all thirteen evaluations failed to reconstruct")


def reconstruct_six(values):
    """Reconstruct T from D_ei and D_ei+ej.

    ``values`` maps the six standard points to their 3x3 D matrices.  In
    characteristic three Delta D(e_i,e_k)=D_i+D_k-D_{i+k} is exactly the
    outer contraction T(e_i,-,e_k), including i=k.
    """
    tensor = [[0] * 27 for _ in range(3)]
    basis_values = [values[p] for p in BASIS_POINTS]
    pair_values = {
        (0, 1): values[(1, 1, 0)],
        (0, 2): values[(1, 0, 1)],
        (1, 2): values[(0, 1, 1)],
    }
    for i in range(3):
        for m in range(3):
            for l in range(3):
                tensor[l][word_idx(i, m, i)] = basis_values[i][m][l]
    for i, k in ((0, 1), (0, 2), (1, 2)):
        for m in range(3):
            for l in range(3):
                value = (basis_values[i][m][l] + basis_values[k][m][l]
                         - pair_values[(i, k)][m][l]) % P
                tensor[l][word_idx(i, m, k)] = value
                tensor[l][word_idx(k, m, i)] = value
    assert validate_tensor(tensor)
    return tensor


def solve_unique(matrix, rhs):
    augmented = [row[:] + [value % P] for row, value in zip(matrix, rhs)]
    rows, pivots = rref(augmented)
    if len([p for p in pivots if p < len(matrix[0])]) != len(matrix[0]):
        raise ValueError("system is not uniquely solvable")
    solution = [0] * len(matrix[0])
    for row, pivot in zip(rows, pivots):
        if pivot < len(solution):
            solution[pivot] = row[-1]
    return solution


def reconstruct_from_points(points, values):
    """Generic exact inversion for any full-rank evaluation family."""
    basis = admissible_basis()
    matrix = evaluation_matrix(points)
    rows = []
    for l in range(3):
        rhs = [values[point][m][l] for point in points for m in range(3)]
        coefficients = solve_unique(matrix, rhs)
        rows.append([
            sum(coefficients[q] * basis[q][w] for q in range(11)) % P
            for w in range(27)
        ])
    assert validate_tensor(rows)
    return rows


def main():
    size, families = minimal_families()
    points = projective_points()
    preferred = tuple(sorted(points.index(p) for p in
                      [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]))
    assert preferred in families
    chosen_indices = preferred
    chosen_points = [points[i] for i in chosen_indices]

    rng = random.Random(SEED)
    for _ in range(250):
        tensor = random_tensor(rng)
        six_values = {point: contract_D(tensor, point)
                      for point in BASIS_POINTS + PAIR_POINTS}
        assert reconstruct_six(six_values) == tensor
        minimal_values = {point: contract_D(tensor, point)
                          for point in chosen_points}
        assert reconstruct_from_points(chosen_points, minimal_values) == tensor

        # The live diagonal is already visible in the basis evaluations:
        # delta[l,i] is the i-th row of D_{e_i}, column l.
        delta = diagonal_matrix(tensor)
        for i, point in enumerate(BASIS_POINTS):
            d = six_values[point]
            assert [d[i][l] for l in range(3)] == [delta[l][i] for l in range(3)]

    result = {
        "prime": P,
        "seed": SEED,
        "random_trials": 250,
        "admissible_dimension_per_relation": 11,
        "tensor_parameter_count": 33,
        "all_projective_evaluation_rank": evaluation_rank(points),
        "six_standard_evaluation_rank": evaluation_rank(BASIS_POINTS + PAIR_POINTS),
        "minimal_number_of_full_matrix_evaluations": size,
        "number_of_minimal_families": len(families),
        "first_minimal_family_indices": list(chosen_indices),
        "first_minimal_family_points": chosen_points,
        "all_minimal_families": [[points[i] for i in combo] for combo in families],
        "diagonal_is_extra_beyond_full_D_values": False,
        "moduli_dimension_heuristic": 16,
        "polarization": "DeltaD(x,z)=D_x+D_z-D_{x+z}=T(x,-,z)",
    }
    out = HERE / "results" / "reconstruction.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k != "all_minimal_families"}, indent=2))


if __name__ == "__main__":
    main()
