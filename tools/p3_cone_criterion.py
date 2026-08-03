#!/usr/bin/env python3
"""Bockstein-cone criterion at p=3 on the verified source tensors.

For every tensor in examples/p3/source-tensors/, this script computes
the Bockstein matrix B (the diagonal coefficients), the reduced linear
cone C_beta = P(ker B), the closed points of the norm-degeneracy scheme
Sigma_D of residue degree at most three, and, at every rank-one cone
point, the transversality map Theta_x.  At each reduced isolated cone
point it constructs the adapted basis of the criterion and checks that
the transformed relation rows have the pivot words XXZ, XYY, XYZ.

Everything is exact linear algebra over F_3, F_9, F_27; the script is
deterministic and writes examples/p3/cone-criterion/report.txt.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "examples" / "p3" / "source-tensors"
OUT = ROOT / "examples" / "p3" / "cone-criterion"

# ----------------------------------------------------------------------
# Finite fields F_{3^d} as coefficient tuples modulo a fixed minimal
# polynomial; the moduli are x^2+1 over F_3 and x^3-x+1 over F_3.

MODULI = {1: (0,), 2: (1, 0), 3: (1, 2, 0)}  # low-degree coefficients


class GF:
    def __init__(self, d):
        self.d = d
        self.mod = MODULI[d]

    def zero(self):
        return (0,) * self.d

    def one(self):
        return (1,) + (0,) * (self.d - 1)

    def from_int(self, n):
        return (n % 3,) + (0,) * (self.d - 1)

    def add(self, a, b):
        return tuple((x + y) % 3 for x, y in zip(a, b))

    def neg(self, a):
        return tuple((-x) % 3 for x in a)

    def mul(self, a, b):
        prod = [0] * (2 * self.d - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    prod[i + j] = (prod[i + j] + x * y) % 3
        for k in range(2 * self.d - 2, self.d - 1, -1):
            c = prod[k]
            if c:
                prod[k] = 0
                for j, m in enumerate(self.mod):
                    prod[k - self.d + j] = (prod[k - self.d + j] - c * m) % 3
        return tuple(prod[: self.d])

    def inv(self, a):
        # brute force is fine in a field of at most 27 elements
        for b in self.elements():
            if self.mul(a, b) == self.one():
                return b
        raise ZeroDivisionError

    def elements(self):
        from itertools import product
        return [tuple(t) for t in product(range(3), repeat=self.d)]

    def frob(self, a):
        return self.mul(self.mul(a, a), a)


def mat_rank(k, rows):
    m = [list(r) for r in rows]
    rank, col, ncols = 0, 0, len(m[0]) if m else 0
    while rank < len(m) and col < ncols:
        piv = next((r for r in range(rank, len(m)) if m[r][col] != k.zero()), None)
        if piv is None:
            col += 1
            continue
        m[rank], m[piv] = m[piv], m[rank]
        ipiv = k.inv(m[rank][col])
        m[rank] = [k.mul(ipiv, v) for v in m[rank]]
        for r in range(len(m)):
            if r != rank and m[r][col] != k.zero():
                c = m[r][col]
                m[r] = [k.add(v, k.neg(k.mul(c, w))) for v, w in zip(m[r], m[rank])]
        rank += 1
        col += 1
    return rank, m


def kernel_basis(k, rows):
    n = len(rows[0])
    rank, rref = mat_rank(k, rows)
    pivots = []
    for r in range(rank):
        pivots.append(next(c for c in range(n) if rref[r][c] != k.zero()))
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for f in free:
        v = [k.zero()] * n
        v[f] = k.one()
        for r, pc in enumerate(pivots):
            v[pc] = k.neg(rref[r][f])
        basis.append(tuple(v))
    return basis


# ----------------------------------------------------------------------
# Tensor contractions.  T[ell][idx] with idx = 9(i-1)+3(j-1)+(k-1) is
# m_{ell,ijk} = <M(chi_i,chi_j,chi_k), e_ell>.

def word_index(i, j, k):
    return 9 * i + 3 * j + k  # zero-based letters


def D_matrix(k, T, x):
    """Matrix of D_x over k: entry (j, ell) = chi_j(D_x(e_ell))."""
    mat = [[k.zero()] * 3 for _ in range(3)]
    for j in range(3):
        for ell in range(3):
            s = k.zero()
            for i in range(3):
                for kk in range(3):
                    c = T[ell][word_index(i, j, kk)] % 3
                    if c:
                        s = k.add(s, k.mul(k.from_int(c), k.mul(x[i], x[kk])))
            mat[j][ell] = s
    return mat


def bockstein_matrix(T):
    return [[T[ell][word_index(i, i, i)] % 3 for i in range(3)] for ell in range(3)]


def cone_value(k, B, x):
    return tuple(
        (lambda s: s)(
            _sum(k, [k.mul(k.from_int(B[ell][i]), x[i]) for i in range(3)])
        )
        for ell in range(3)
    )


def _sum(k, vals):
    s = k.zero()
    for v in vals:
        s = k.add(s, v)
    return s


def proj_points(k):
    pts, seen = [], set()
    for v in _all_vectors(k):
        if all(c == k.zero() for c in v):
            continue
        n = _normalize(k, v)
        if n not in seen:
            seen.add(n)
            pts.append(n)
    return pts


def _all_vectors(k):
    from itertools import product
    els = k.elements()
    return [tuple(t) for t in product(els, repeat=3)]


def _normalize(k, v):
    lead = next(c for c in v if c != k.zero())
    il = k.inv(lead)
    return tuple(k.mul(il, c) for c in v)


def frob_point(k, v):
    return _normalize(k, tuple(k.frob(c) for c in v))


def theta_rank(k, T, x, Dx):
    """Rank of Theta_x at a rank-one cone point (adapted-datum route)."""
    cols = [tuple(Dx[j][ell] for j in range(3)) for ell in range(3)]
    imrank, _ = mat_rank(k, [list(c) for c in cols])
    assert imrank == 1
    # kernel of D_x (inputs e with D_x e = 0)
    ker = kernel_basis(k, [[Dx[j][ell] for ell in range(3)] for j in range(3)])
    assert len(ker) == 2
    # image line and a functional y with y(im)=0, y not proportional to x
    im = next(c for c in cols if any(v != k.zero() for v in c))
    ycands = kernel_basis(k, [list(im)])
    y = next(
        c for c in ycands
        if mat_rank(k, [list(x), list(c)])[0] == 2
    )
    # representatives v1, v2 of a basis of V/kx
    reps = []
    for v in _std_vectors(k):
        if mat_rank(k, [list(x)] + [list(r) for r in reps] + [list(v)])[0] == 2 + len(reps):
            reps.append(v)
        if len(reps) == 2:
            break
    rows = []
    for v in reps:
        xv = tuple(k.add(a, b) for a, b in zip(x, v))
        Dv, Dxv = D_matrix(k, T, v), D_matrix(k, T, xv)
        delta = [[k.add(Dx[j][l], k.add(Dv[j][l], k.neg(Dxv[j][l]))) for l in range(3)]
                 for j in range(3)]
        row = []
        for e in ker:
            val = k.zero()
            for j in range(3):
                comp = _sum(k, [k.mul(delta[j][l], e[l]) for l in range(3)])
                val = k.add(val, k.mul(y[j], comp))
            row.append(k.neg(val))
        rows.append(row)
    r, _ = mat_rank(k, rows)
    return r


def _std_vectors(k):
    outs = []
    for i in range(3):
        v = [k.zero()] * 3
        v[i] = k.one()
        outs.append(tuple(v))
    return outs


DEGLEX = sorted(range(27), key=lambda n: (-(n // 9 == 0) * 0, ), reverse=False)


def deglex_order():
    """Word indices sorted descending for X>Y>Z (letters 0>1>2)."""
    words = [(i, j, kk) for i in range(3) for j in range(3) for kk in range(3)]
    words.sort()  # lexicographic with 0 (=X) smallest tuple value = largest letter
    return [word_index(*w) for w in words]


def anick_heads(k, T, basis):
    """Pivot words of the transformed relation rows for character basis
    rows A (each an element of V tensor k)."""
    A = basis
    newT = []
    for ell in range(3):
        row = [k.zero()] * 27
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    s = k.zero()
                    for i in range(3):
                        for j in range(3):
                            for kk in range(3):
                                coef = T[ell][word_index(i, j, kk)] % 3
                                if coef:
                                    s = k.add(
                                        s,
                                        k.mul(
                                            k.from_int(coef),
                                            k.mul(A[a][i], k.mul(A[b][j], A[c][kk])),
                                        ),
                                    )
                    row[word_index(a, b, c)] = s
        newT.append(row)
    order = deglex_order()
    perm = [[r[c] for c in order] for r in newT]
    rank, rref = mat_rank(k, perm)
    heads = []
    for r in range(rank):
        c = next(c for c in range(27) if rref[r][c] != k.zero())
        heads.append(order[c])
    return heads, rank


def word_name(idx):
    letters = "XYZ"
    i, r = divmod(idx, 9)
    j, kk = divmod(r, 3)
    return letters[i] + letters[j] + letters[kk]


def adapted_basis(k, x, Dx):
    cols = [tuple(Dx[j][ell] for j in range(3)) for ell in range(3)]
    c3 = next(c for c in cols if any(v != k.zero() for v in c))
    xperp = kernel_basis(k, [list(x)])
    c2 = next(c for c in xperp if mat_rank(k, [list(c3), list(c)])[0] == 2)
    c1 = next(v for v in _std_vectors(k)
              if mat_rank(k, [list(c1v) for c1v in (c3, c2, v)])[0] == 3)
    # dual basis (x', y, z) of (c1, c2, c3): rows of inverse of matrix with
    # columns c1, c2, c3
    M = [[c1[j], c2[j], c3[j]] for j in range(3)]
    inv = _mat_inverse(k, M)
    return [tuple(inv[r]) for r in range(3)]


def _mat_inverse(k, M):
    n = len(M)
    aug = [list(M[r]) + [k.one() if c == r else k.zero() for c in range(n)]
           for r in range(n)]
    rank, rref = mat_rank(k, aug)
    assert rank == n
    return [row[n:] for row in rref]


def main():
    OUT.mkdir(exist_ok=True)
    fields = GF(1), GF(2), GF(3)
    lines = []
    for tdir in sorted(SRC.iterdir(), key=lambda p: int(p.name.split("-")[1])):
        data = json.loads((tdir / "tensor.json").read_text())
        T = data["tensor_3_by_27"]
        disc = data["discriminant"]
        B = bockstein_matrix(T)
        k1 = fields[0]
        brank, _ = mat_rank(k1, [[k1.from_int(B[r][c]) for c in range(3)]
                                 for r in range(3)])
        lines.append(f"D_K = {disc}")
        lines.append(f"  Bockstein rank {brank}; cone dimension {2 - brank}"
                     if brank < 3 else "  Bockstein rank 3; cone empty")
        found = []
        seen_orbits = set()
        for k in fields:
            for pt in proj_points(k):
                Dx = D_matrix(k, T, pt)
                rank, _ = mat_rank(k, [[Dx[j][l] for l in range(3)] for j in range(3)])
                if rank > 1:
                    continue
                orbit = []
                q = pt
                while q not in orbit:
                    orbit.append(q)
                    q = frob_point(k, q)
                if len(orbit) != k.d:
                    continue  # counted in the smaller field
                key = (k.d, tuple(sorted(orbit)))
                if key in seen_orbits:
                    continue
                seen_orbits.add(key)
                on_cone = all(v == k.zero() for v in cone_value(k, B, pt))
                entry = {
                    "degree": k.d, "point": pt, "rank": rank, "cone": on_cone,
                }
                if rank == 1 and on_cone:
                    tr = theta_rank(k, T, pt, Dx)
                    entry["theta_rank"] = tr
                    entry["reduced_isolated"] = tr == 2
                    if tr == 2:
                        basis = adapted_basis(k, pt, Dx)
                        heads, rk = anick_heads(k, T, basis)
                        entry["heads"] = sorted(word_name(h) for h in heads)
                        entry["row_rank"] = rk
                found.append(entry)
        if not found:
            lines.append("  Sigma_D has no closed point of degree <= 3 "
                         "with rank at most one")
        for e in found:
            desc = (f"  point degree {e['degree']}, rank {e['rank']}, "
                    f"on cone: {e['cone']}")
            if "theta_rank" in e:
                desc += (f", rk Theta = {e['theta_rank']}, reduced isolated: "
                         f"{e['reduced_isolated']}")
            if "heads" in e:
                desc += (f"; adapted pivot words {'/'.join(e['heads'])} "
                         f"(row rank {e['row_rank']})")
            lines.append(desc)
        lines.append("")
    report = "\n".join(lines) + "\n"
    (OUT / "report.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
