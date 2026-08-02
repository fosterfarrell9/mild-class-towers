#!/usr/bin/env python3
"""Finite algebra for cubic Zassenhaus relations in characteristic three.

The word order is X_i X_j X_k with k fastest.  A relation row ``t`` is
admissible when

    t[i,j,k] = t[k,j,i],
    t[i,j,k] + t[j,k,i] + t[k,i,j] = 0.

In characteristic three the second identity imposes no condition on
``t[i,i,i]``.  Thus the admissible space is the degree-three piece of the
free restricted Lie algebra: ordinary L_3 plus x^3,y^3,z^3.

Everything in this module is deterministic finite linear algebra; it does
not import or invoke number-field functionality.
"""

from __future__ import annotations

from itertools import product
import random

P = 3
N = 3
WORDS = N ** 3
LETTERS = "xyz"


def word_idx(i: int, j: int, k: int) -> int:
    return (i * N + j) * N + k


def decode_word(index: int) -> tuple[int, int, int]:
    return index // 9, (index // 3) % 3, index % 3


def word_name(index: int) -> str:
    return "".join(LETTERS[i] for i in decode_word(index))


def rref(matrix: list[list[int]], column_order=None):
    """RREF over F_3, returning nonzero rows and pivot columns.

    ``column_order`` may specify the order in which pivots are sought.  This
    is used for leading monomials, where largest words are visited first.
    """
    m = [[x % P for x in row] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if rows else 0
    order = list(range(cols)) if column_order is None else list(column_order)
    pivots: list[int] = []
    r = 0
    for c in order:
        pivot = next((i for i in range(r, rows) if m[i][c]), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        inv = pow(m[r][c], P - 2, P)
        m[r] = [(v * inv) % P for v in m[r]]
        for i in range(rows):
            if i != r and m[i][c]:
                f = m[i][c]
                m[i] = [(a - f * b) % P for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return m[:r], pivots


def rank(matrix: list[list[int]]) -> int:
    return len(rref(matrix)[0])


def kernel(matrix: list[list[int]]) -> list[list[int]]:
    """Basis of the right kernel, returned as a list of column-vectors."""
    rows, pivots = rref(matrix)
    n = len(matrix[0]) if matrix else 0
    free = [c for c in range(n) if c not in pivots]
    out = []
    for fc in free:
        v = [0] * n
        v[fc] = 1
        for rr, pc in enumerate(pivots):
            v[pc] = (-rows[rr][fc]) % P
        out.append(v)
    return out


def canonical_rowspace(matrix: list[list[int]]) -> tuple[tuple[int, ...], ...]:
    rows, _ = rref(matrix)
    return tuple(tuple(row) for row in rows)


def admissibility_constraints() -> list[list[int]]:
    constraints: list[list[int]] = []
    for i, j, k in product(range(N), repeat=3):
        row = [0] * WORDS
        row[word_idx(i, j, k)] += 1
        row[word_idx(k, j, i)] -= 1
        constraints.append([v % P for v in row])

        row = [0] * WORDS
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            row[word_idx(a, b, c)] += 1
        constraints.append([v % P for v in row])
    return constraints


def nullspace_basis(constraints: list[list[int]]) -> list[list[int]]:
    rows, pivots = rref(constraints)
    free = [c for c in range(WORDS) if c not in pivots]
    basis = []
    for fc in free:
        v = [0] * WORDS
        v[fc] = 1
        for rr, pc in enumerate(pivots):
            v[pc] = (-rows[rr][fc]) % P
        basis.append(v)
    return basis


def admissible_basis() -> list[list[int]]:
    basis = nullspace_basis(admissibility_constraints())
    assert len(basis) == 11
    return basis


def bracket3(a: int, b: int, c: int) -> list[int]:
    """Associative expansion of [a,[b,c]]."""
    v = [0] * WORDS
    for i, j, k, coefficient in (
        (a, b, c, 1), (a, c, b, -1),
        (b, c, a, -1), (c, b, a, 1),
    ):
        q = word_idx(i, j, k)
        v[q] = (v[q] + coefficient) % P
    return v


def ordinary_lie_basis() -> list[list[int]]:
    """A Hall-type basis of the ordinary free Lie piece L_3."""
    triples = [
        (0, 0, 1), (0, 0, 2), (1, 0, 1), (1, 1, 2),
        (2, 0, 2), (2, 1, 2), (0, 1, 2), (1, 0, 2),
    ]
    rows = [bracket3(*triple) for triple in triples]
    assert rank(rows) == 8
    return rows


def restricted_lie_basis() -> list[list[int]]:
    basis = ordinary_lie_basis()
    for i in range(N):
        cube = [0] * WORDS
        cube[word_idx(i, i, i)] = 1
        basis.append(cube)
    assert rank(basis) == 11
    return basis


def validate_row(row: list[int]) -> bool:
    if len(row) != WORDS:
        return False
    for i, j, k in product(range(N), repeat=3):
        if (row[word_idx(i, j, k)] - row[word_idx(k, j, i)]) % P:
            return False
        if sum(row[word_idx(*q)] for q in
               ((i, j, k), (j, k, i), (k, i, j))) % P:
            return False
    return True


def validate_tensor(tensor: list[list[int]], require_rank_three=True) -> bool:
    return (len(tensor) == 3
            and all(validate_row(row) for row in tensor)
            and (not require_rank_three or rank(tensor) == 3))


def projective_points() -> list[tuple[int, int, int]]:
    """The thirteen normalized points of P^2(F_3)."""
    return ([(1, b, c) for b in range(P) for c in range(P)]
            + [(0, 1, c) for c in range(P)]
            + [(0, 0, 1)])


def contract_row(row: list[int], point) -> list[int]:
    """D_x for one relation row, as a vector in C=V^*.

    At p=3 the usual contraction is already D_x because -2=1.
    """
    return [
        sum(point[i] * point[k] * row[word_idx(i, m, k)]
            for i in range(N) for k in range(N)) % P
        for m in range(N)
    ]


def contract_D(tensor: list[list[int]], point) -> list[list[int]]:
    """Return the 3x3 matrix D_x, rows in C and columns in E."""
    columns = [contract_row(row, point) for row in tensor]
    return [[columns[l][m] for l in range(3)] for m in range(3)]


def diagonal_matrix(tensor: list[list[int]]) -> list[list[int]]:
    """Matrix delta with delta[l,i]=T[l,i,i,i]."""
    return [[row[word_idx(i, i, i)] for i in range(N)] for row in tensor]


def euler_defect(tensor: list[list[int]], point) -> list[int]:
    """x^t D_x; this equals delta(x^[3]), not zero in general."""
    d = contract_D(tensor, point)
    return [sum(point[m] * d[m][l] for m in range(N)) % P
            for l in range(3)]


def frobenius_diagonal(tensor: list[list[int]], point) -> list[int]:
    delta = diagonal_matrix(tensor)
    return [sum(delta[l][i] * point[i] ** 3 for i in range(N)) % P
            for l in range(3)]


def determinant3(matrix: list[list[int]]) -> int:
    a = matrix
    return (a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
            - a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
            + a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0])) % P


def rank_signature(tensor: list[list[int]]) -> tuple[int, int, int, int]:
    ranks = [rank(contract_D(tensor, point)) for point in projective_points()]
    return tuple(ranks.count(r) for r in range(4))


def random_tensor(rng: random.Random) -> list[list[int]]:
    basis = admissible_basis()
    while True:
        coeffs = [[rng.randrange(P) for _ in basis] for _ in range(3)]
        tensor = [[sum(c*basis[q][w] for q, c in enumerate(row)) % P
                   for w in range(WORDS)] for row in coeffs]
        if rank(tensor) == 3:
            assert validate_tensor(tensor)
            return tensor


def hall_and_cube_basis():
    lie = ordinary_lie_basis()
    names = [
        "[x,[x,y]]", "[x,[x,z]]", "[y,[x,y]]", "[y,[y,z]]",
        "[z,[x,z]]", "[z,[y,z]]", "[x,[y,z]]", "[y,[x,z]]",
    ]
    named = list(zip(names, lie))
    for i, name in enumerate(LETTERS):
        cube = [0] * WORDS
        cube[word_idx(i, i, i)] = 1
        named.append((f"{name}^3", cube))
    return named


def gl3_matrices():
    """Yield all 11232 elements of GL_3(F_3), deterministically."""
    vectors = [v for v in product(range(P), repeat=3) if any(v)]
    vectors.sort(key=lambda v: (sum(x != 0 for x in v), v))
    for a in vectors:
        for b in vectors:
            for c in vectors:
                matrix = [[a[i], b[i], c[i]] for i in range(3)]
                if determinant3(matrix):
                    yield matrix


def transform_row(row: list[int], matrix: list[list[int]]) -> list[int]:
    """Substitute X_i -> sum_a matrix[a][i] X_a."""
    first = [[[sum(matrix[a][i] * row[word_idx(i, j, k)]
                    for i in range(3)) % P
               for k in range(3)] for j in range(3)] for a in range(3)]
    second = [[[sum(matrix[b][j] * first[a][j][k]
                     for j in range(3)) % P
                for k in range(3)] for b in range(3)] for a in range(3)]
    return [sum(matrix[c][k] * second[a][b][k] for k in range(3)) % P
            for a, b, c in product(range(3), repeat=3)]


def transform_tensor(tensor, matrix):
    return [transform_row(row, matrix) for row in tensor]


def combinatorially_free(words: list[int]) -> bool:
    decoded = [decode_word(word) for word in words]
    if len(set(words)) != len(words):
        return False
    for left in decoded:
        for right in decoded:
            for overlap in (1, 2):
                if left[:overlap] == right[3-overlap:]:
                    return False
    return True


def leading_words(tensor: list[list[int]]) -> list[int]:
    """Leading words after relation-basis reduction for x<y<z lex."""
    _, pivots = rref(tensor, range(WORDS - 1, -1, -1))
    assert len(pivots) == 3
    return pivots


def find_anick_witness(tensor: list[list[int]]):
    """Exhaust GL_3(F_3) for a rational Anick witness.

    Testing one fixed lexicographic letter order is exhaustive: changing the
    letter order is a permutation matrix and is already absorbed into GL_3.
    """
    tested = 0
    for matrix in gl3_matrices():
        tested += 1
        transformed = transform_tensor(tensor, matrix)
        leaders = leading_words(transformed)
        if combinatorially_free(leaders):
            return {
                "found": True,
                "tested": tested,
                "matrix": matrix,
                "leaders": [word_name(q) for q in leaders],
            }
    assert tested == 11232
    return {"found": False, "tested": tested, "matrix": None, "leaders": []}


def self_check() -> dict:
    identities = admissible_basis()
    restricted = restricted_lie_basis()
    assert canonical_rowspace(identities) == canonical_rowspace(restricted)
    assert len(projective_points()) == 13
    assert sum(1 for _ in gl3_matrices()) == 11232
    rng = random.Random(20260801)
    for _ in range(25):
        tensor = random_tensor(rng)
        for point in projective_points():
            assert euler_defect(tensor, point) == frobenius_diagonal(tensor, point)
    return {
        "identity_dimension": len(identities),
        "ordinary_lie_dimension": len(ordinary_lie_basis()),
        "restricted_lie_dimension": len(restricted),
        "projective_points": 13,
        "gl3_order": 11232,
        "seed": 20260801,
    }


if __name__ == "__main__":
    print(self_check())
