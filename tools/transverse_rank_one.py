#!/usr/bin/env python3
"""Transverse rank-one certificates for the computed fields.

For each field the quadratic secondary-norm family D(u,v,w) is rebuilt by
polarization from the six verified matrices stored in the repository (the
committed result.gp records; for the principal example the file
certificate/K-2800905-p5/secondary-norms.gp exported from its arithmetic
certificate by tools/export_secondary_norms.gp).  The transverse rank-one
criterion of the paper,

    rank(D_x) = 1   and   det(B_x) != 0,
    (B_x)_{ij} = y(DeltaD(x, v_i)(e_j)),
    e_2,e_3 basis of ker D_x,  v_1,v_2 basis of V_k/kx,
    y in (im D_x)-perp, y not in k x,
    DeltaD(x,v) = D_x + D_v - D_{x+v},

is tested at all 31 rational points of P^2(F_5) and at every non-rational
closed point of the rank-drop scheme Sigma_D = {rank D_x <= 1}.  Closed
points are located exactly: per affine chart the radical of the ideal of
2x2 minors is computed with Singular (lex order), the univariate
eliminant's roots over F_{5^d}, d <= 6, are found by exhaustive scanning,
and every certificate is direct linear algebra over an explicitly
presented field k = F_5[t]/(f(t)).

Internal verification before any certificate is produced:
  * the family reproduces the six stored matrices of every field;
  * the Euler identity x . D_x = 0 holds symbolically (so the family is
    a bundle map E (x) O -> Omega^1(3), and rank D_x <= 2 throughout);
  * for the eight fields whose result record stores the 3 x 27 cubic
    relation matrix, the family reconstructs it exactly (convention
    lock across all fields);
  * at every certified point a Jacobian tangent-space computation must
    agree with det(B_x) != 0;
  * PASS/FAIL is checked to be invariant under seeded random changes of
    the kernel and quotient bases.

Deterministic.  Requires Singular.  Usage, from the repository root:
    python3 tools/transverse_rank_one.py
Writes examples/p5/transverse-rank-one/{certificates.json,report.txt}.
"""

from __future__ import annotations

import json
import random
import re
import subprocess
import sys
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "examples" / "p5" / "transverse-rank-one"

P = 5
SEED = 20260728
SIX = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]

FIELDS = [
    (-11203620, ROOT / "certificate/K-2800905-p5/secondary-norms.gp"),
    (-18397407, ROOT / "examples/p5/D-18397407/result.gp"),
    (-27960639, ROOT / "examples/p5/batch-block0-01/D-27960639/result.gp"),
    (-35663739, ROOT / "examples/p5/batch-block0-01/D-35663739/result.gp"),
    (-51213139, ROOT / "examples/p5/batch-block0-01/D-51213139/result.gp"),
    (-54319112, ROOT / "examples/p5/batch-block0-01/D-54319112/result.gp"),
    (-61040707, ROOT / "examples/p5/batch-block0-01/D-61040707/result.gp"),
    (-65818135, ROOT / "examples/p5/batch-block0-01/D-65818135/result.gp"),
    (-75949255, ROOT / "examples/p5/batch-block0-01/D-75949255/result.gp"),
    (-145367147, ROOT / "examples/p5/batch-block0-01/D-145367147/result.gp"),
    (-207666763, ROOT / "examples/p5/batch-block0-01/D-207666763/result.gp"),
    (-104545864, ROOT / "examples/p5/further/D-104545864/result.gp"),
    (-123779560, ROOT / "examples/p5/further/D-123779560/result.gp"),
    (-126740891, ROOT / "examples/p5/further/D-126740891/result.gp"),
    (-89017304, ROOT / "examples/p5/further/D-89017304/result.gp"),
    (-89218664, ROOT / "examples/p5/further/D-89218664/result.gp"),
    (-90903207, ROOT / "examples/p5/further/D-90903207/result.gp"),
    (-93121640, ROOT / "examples/p5/further/D-93121640/result.gp"),
    (-123482119, ROOT / "examples/p5/further/D-123482119/result.gp"),
    (-106660295, ROOT / "examples/p5/further/D-106660295/result.gp"),
    (-77778287, ROOT / "examples/p5/further/D-77778287/result.gp"),
]


# ------------------------------------------------------------ input parsing

def parse_gp_matrices(text, key):
    match = re.search(r'"%s",\s*\[(.*?)\]\]' % key, text, re.S)
    if not match:
        raise SystemExit(f"missing {key}")
    body = match.group(1) + "]"
    mats = []
    for chunk in re.findall(r"\[([^\[\]]*)\]", body):
        mats.append([[int(v) % P for v in row.split(",")]
                     for row in chunk.split(";")])
    return mats


def parse_relation_matrix(text):
    match = re.search(r'"cubic_relation_matrix",\s*\[(.*?)\]\]', text, re.S)
    if not match:
        return None
    rows = [[int(v) % P for v in row.split(",")]
            for row in match.group(1).split(";")]
    assert len(rows) == 3 and all(len(r) == 27 for r in rows)
    return rows


# --------------------------------------------------- F_5 linear algebra

def rref(matrix):
    m = [row[:] for row in matrix]
    rows = len(m)
    cols = len(m[0]) if rows else 0
    pivots = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if m[i][c] % P), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        inv = pow(m[r][c], P - 2, P)
        m[r] = [(v * inv) % P for v in m[r]]
        for i in range(rows):
            if i != r and m[i][c] % P:
                f = m[i][c]
                m[i] = [(m[i][k] - f * m[r][k]) % P for k in range(cols)]
        pivots.append(c)
        r += 1
    return [row for row in m[:r]], pivots


def int_rank(matrix):
    return len(rref(matrix)[0])


# --------------------------------------------- polarization and conventions

def mat_add_int(a, b):
    return [[(x + y) % P for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def mat_scale_int(a, s):
    return [[(x * s) % P for x in row] for row in a]


PAIR_KEY = {(0, 1): 0, (0, 2): 1, (1, 2): 2}


def D_of(basis, pairs, lam):
    result = [[0] * 3 for _ in range(3)]
    for i in range(3):
        c = (lam[i] * lam[i]) % P
        if c:
            result = mat_add_int(result, mat_scale_int(basis[i], c))
    for (i, j), pk in PAIR_KEY.items():
        c = (lam[i] * lam[j]) % P
        if c:
            cross = mat_add_int(
                pairs[pk],
                mat_scale_int(mat_add_int(basis[i], basis[j]), P - 1))
            result = mat_add_int(result, mat_scale_int(cross, c))
    return result


def delta_D_int(basis, pairs, i, k):
    if i == k:
        return mat_scale_int(basis[i], P - 2)
    cross = mat_add_int(basis[i], basis[k])
    return mat_add_int(
        cross, mat_scale_int(pairs[PAIR_KEY[(min(i, k), max(i, k))]],
                             P - 1))


def word_matrix(basis, pairs):
    T = [[0] * 27 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                dd = delta_D_int(basis, pairs, i, k)
                for l in range(3):
                    T[l][(i * 3 + j) * 3 + k] = dd[j][l]
    return T


# -------------------------------------------- symbolic family and minors

MONO_SQ = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
MONO_CR = {(0, 1): (1, 1, 0), (0, 2): (1, 0, 1), (1, 2): (0, 1, 1)}


def poly_add(a, b):
    out = dict(a)
    for mono, c in b.items():
        out[mono] = (out.get(mono, 0) + c) % P
        if not out[mono]:
            del out[mono]
    return out


def poly_scale(a, s):
    return {m: (c * s) % P for m, c in a.items() if (c * s) % P}


def poly_mul(a, b):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            out[m] = (out.get(m, 0) + c1 * c2) % P
    return {m: c for m, c in out.items() if c}


def symbolic_D(basis, pairs):
    D = [[{} for _ in range(3)] for _ in range(3)]
    for m in range(3):
        for l in range(3):
            entry = {}
            for i in range(3):
                if basis[i][m][l]:
                    entry = poly_add(entry, {MONO_SQ[i]: basis[i][m][l]})
            for (i, j), mono in MONO_CR.items():
                c = (pairs[PAIR_KEY[(i, j)]][m][l] - basis[i][m][l]
                     - basis[j][m][l]) % P
                if c:
                    entry = poly_add(entry, {mono: c})
            D[m][l] = entry
    return D


def euler_check(D):
    lam = [{(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}]
    for l in range(3):
        total = {}
        for m in range(3):
            total = poly_add(total, poly_mul(lam[m], D[m][l]))
        if total:
            return False
    return True


def minors2(D):
    from itertools import combinations
    out = []
    for r1, r2 in combinations(range(3), 2):
        for c1, c2 in combinations(range(3), 2):
            out.append(poly_add(
                poly_mul(D[r1][c1], D[r2][c2]),
                poly_scale(poly_mul(D[r1][c2], D[r2][c1]), P - 1)))
    return out


# ------------------------------------------------------------ F_5[t] basics

def f5_trim(a):
    while a and a[-1] % P == 0:
        a.pop()
    return a


def f5_mul(a, b):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x % P:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % P
    return f5_trim(out)


def f5_rem(a, m):
    a = [x % P for x in a]
    dm = len(m) - 1
    inv = pow(m[-1], P - 2, P)
    while f5_trim(a) and len(a) - 1 >= dm:
        c = (a[-1] * inv) % P
        shift = len(a) - 1 - dm
        for i, y in enumerate(m):
            a[shift + i] = (a[shift + i] - c * y) % P
        f5_trim(a)
    return a


def f5_divmod(a, b):
    a = [x % P for x in a]
    b = f5_trim([x % P for x in b])
    inv = pow(b[-1], P - 2, P)
    q = [0] * max(1, len(a) - len(b) + 1)
    while f5_trim(a) and len(a) - 1 >= len(b) - 1:
        c = (a[-1] * inv) % P
        shift = len(a) - len(b)
        q[shift] = c
        for i, y in enumerate(b):
            a[shift + i] = (a[shift + i] - c * y) % P
    return f5_trim(q), f5_trim(a)


def f5_gcd(a, b):
    a, b = list(a), list(b)
    while f5_trim(b):
        a, b = b, f5_rem(a, b)
    return f5_trim(a)


def zip_pad(a, b):
    n = max(len(a), len(b))
    return zip(a + [0] * (n - len(a)), b + [0] * (n - len(b)))


def f5_powmod_x(q, m):
    result = [1]
    base = [0, 1]
    while q:
        if q & 1:
            result = f5_rem(f5_mul(result, base), m)
        base = f5_rem(f5_mul(base, base), m)
        q >>= 1
    return result


def is_irreducible(m):
    d = len(m) - 1
    if f5_trim([(x - y) % P for x, y in
                zip_pad(f5_powmod_x(P ** d, m), [0, 1])]):
        return False
    for pr in {q for q in (2, 3, 5) if d % q == 0}:
        g = f5_gcd([(x - y) % P for x, y in
                    zip_pad(f5_powmod_x(P ** (d // pr), m), [0, 1])], m)
        if len(g) - 1 > 0:
            return False
    return True


def minimal_irreducible(d):
    if d == 1:
        return [0, 1]
    for coeffs in product(range(P), repeat=d):
        m = list(coeffs) + [1]
        if is_irreducible(m):
            return m
    raise AssertionError


# ------------------------------------------------------------------ GF(5^d)

class GF:
    def __init__(self, d):
        self.d = d
        self.modulus = minimal_irreducible(d)
        self.zero = (0,) * d
        self.one = tuple([1] + [0] * (d - 1))

    def from_int(self, c):
        return tuple([c % P] + [0] * (self.d - 1))

    def elem(self, coeffs):
        c = list(coeffs) + [0] * (self.d - len(coeffs))
        return tuple(x % P for x in c[: self.d])

    def add(self, a, b):
        return tuple((x + y) % P for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple((x - y) % P for x, y in zip(a, b))

    def neg(self, a):
        return tuple((-x) % P for x in a)

    def mul(self, a, b):
        return self.elem(f5_rem(f5_mul(list(a), list(b)), self.modulus))

    def inv(self, a):
        r0, r1 = self.modulus[:], f5_trim(list(a))
        if not r1:
            raise ZeroDivisionError
        s0, s1 = [], [1]
        while True:
            f5_trim(r1)
            if len(r1) - 1 <= 0:
                break
            q, r = f5_divmod(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, f5_trim([(x - y) % P for x, y in
                                  zip_pad(s0, f5_mul(q, s1))])
        c = pow(r1[0], P - 2, P)
        return self.elem([x * c % P for x in s1])

    def frob(self, a):
        out = a
        for _ in range(P - 1):
            out = self.mul(out, a)
        return out

    def is_zero(self, a):
        return not any(a)

    def elements(self):
        for coeffs in product(range(P), repeat=self.d):
            yield tuple(coeffs)

    def poly_eval(self, coeffs_f5, a):
        acc = self.zero
        for c in reversed(coeffs_f5):
            acc = self.add(self.mul(acc, a), self.from_int(c))
        return acc


# --------------------------------------------- linear algebra over GF

def mat_rank_ker(field, M):
    rows = [list(r) for r in M]
    nr, nc = len(rows), len(rows[0])
    pivots = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr)
                    if not field.is_zero(rows[i][c])), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = field.inv(rows[r][c])
        rows[r] = [field.mul(x, inv) for x in rows[r]]
        for i in range(nr):
            if i != r and not field.is_zero(rows[i][c]):
                f = rows[i][c]
                rows[i] = [field.sub(x, field.mul(f, y))
                           for x, y in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(nc) if c not in pivots]
    kernel = []
    for fc in free:
        v = [field.zero] * nc
        v[fc] = field.one
        for i, pc in enumerate(pivots):
            v[pc] = field.neg(rows[i][fc])
        kernel.append(v)
    return r, kernel


def reduce_sum(field, xs):
    acc = field.zero
    for x in xs:
        acc = field.add(acc, x)
    return acc


def mat_apply(field, M, v):
    return [reduce_sum(field, [field.mul(M[i][j], v[j])
                               for j in range(len(v))])
            for i in range(len(M))]


def dot(field, a, b):
    return reduce_sum(field, [field.mul(x, y) for x, y in zip(a, b)])


def proportional(field, a, b):
    piv = next(i for i, x in enumerate(b) if not field.is_zero(x))
    if field.is_zero(a[piv]):
        return all(field.is_zero(x) for x in a)
    c = field.mul(a[piv], field.inv(b[piv]))
    return all(field.is_zero(field.sub(x, field.mul(c, y)))
               for x, y in zip(a, b))


# -------------------------------------------------- D-family over GF

def D_matrix(field, basis, pairs, lam):
    def scal(c_int, s):
        return field.mul(field.from_int(c_int), s)
    M = [[field.zero] * 3 for _ in range(3)]
    for i in range(3):
        li2 = field.mul(lam[i], lam[i])
        for m in range(3):
            for l in range(3):
                M[m][l] = field.add(M[m][l], scal(basis[i][m][l], li2))
    for (i, j), pk in PAIR_KEY.items():
        lij = field.mul(lam[i], lam[j])
        for m in range(3):
            for l in range(3):
                c = (pairs[pk][m][l] - basis[i][m][l]
                     - basis[j][m][l]) % P
                M[m][l] = field.add(M[m][l], scal(c, lij))
    return M


def delta_D(field, basis, pairs, x, v):
    Dx = D_matrix(field, basis, pairs, x)
    Dv = D_matrix(field, basis, pairs, v)
    Dxv = D_matrix(field, basis, pairs,
                   [field.add(a, b) for a, b in zip(x, v)])
    return [[field.sub(field.add(Dx[m][l], Dv[m][l]), Dxv[m][l])
             for l in range(3)] for m in range(3)]


def certificate_at(field, basis, pairs, x, rng=None):
    Dx = D_matrix(field, basis, pairs, x)
    rank, ker = mat_rank_ker(field, Dx)
    if rank != 1:
        return {"rank": rank}
    cols = [[Dx[m][l] for m in range(3)] for l in range(3)]
    img = next(c for c in cols if any(not field.is_zero(t) for t in c))
    _, ann = mat_rank_ker(field, [img])
    y = next(v for v in ann if not proportional(field, v, x))
    pivot = next(i for i, t in enumerate(x) if not field.is_zero(t))
    vs = [[field.one if i == m else field.zero for i in range(3)]
          for m in range(3) if m != pivot]
    e2, e3 = ker
    B = [[dot(field, y, mat_apply(field, delta_D(field, basis, pairs,
                                                 list(x), v), e))
          for e in (e2, e3)] for v in vs]
    detB = field.sub(field.mul(B[0][0], B[1][1]),
                     field.mul(B[0][1], B[1][0]))
    result = {
        "rank": 1, "kernel": [e2, e3], "image": img, "vs": vs, "y": y,
        "B": B, "detB": detB, "transverse": not field.is_zero(detB),
    }
    if rng is not None:
        for _ in range(3):
            g = [[field.from_int(rng.randrange(P)) for _ in range(2)]
                 for _ in range(2)]
            if field.is_zero(field.sub(field.mul(g[0][0], g[1][1]),
                                       field.mul(g[0][1], g[1][0]))):
                continue
            e2b = [field.add(field.mul(g[0][0], a), field.mul(g[0][1], b))
                   for a, b in zip(e2, e3)]
            e3b = [field.add(field.mul(g[1][0], a), field.mul(g[1][1], b))
                   for a, b in zip(e2, e3)]
            Bb = [[dot(field, y, mat_apply(
                field, delta_D(field, basis, pairs, list(x), v), e))
                for e in (e2b, e3b)] for v in vs]
            db = field.sub(field.mul(Bb[0][0], Bb[1][1]),
                           field.mul(Bb[0][1], Bb[1][0]))
            assert field.is_zero(db) == field.is_zero(detB), \
                "basis-invariance violated"
    return result


# ------------------------------------------- Jacobian tangent cross-check

def substitute_chart(mono_poly, chart):
    keep = [i for i in range(3) if i != chart]
    out = {}
    for mono, c in mono_poly.items():
        key = (mono[keep[0]], mono[keep[1]])
        out[key] = (out.get(key, 0) + c) % P
    return {k: v for k, v in out.items() if v}


def partial(poly2, var):
    out = {}
    for (a, b), c in poly2.items():
        e = (a, b)[var]
        if e:
            key = (a - 1, b) if var == 0 else (a, b - 1)
            out[key] = (out.get(key, 0) + c * e) % P
    return {k: v for k, v in out.items() if v}


def eval2(field, poly2, s, t):
    acc = field.zero
    for (a, b), c in poly2.items():
        term = field.from_int(c)
        for _ in range(a):
            term = field.mul(term, s)
        for _ in range(b):
            term = field.mul(term, t)
        acc = field.add(acc, term)
    return acc


def tangent_dim(field, minors, x):
    pivot = next(i for i, t in enumerate(x) if not field.is_zero(t))
    inv = field.inv(x[pivot])
    xn = [field.mul(t, inv) for t in x]
    keep = [i for i in range(3) if i != pivot]
    s, t = xn[keep[0]], xn[keep[1]]
    rows = []
    for m in minors:
        m2 = substitute_chart(m, pivot)
        assert field.is_zero(eval2(field, m2, s, t)), "point not on scheme"
        rows.append([eval2(field, partial(m2, 0), s, t),
                     eval2(field, partial(m2, 1), s, t)])
    rank, _ = mat_rank_ker(field, rows)
    return 2 - rank


# ------------------------------------------- locating non-rational points

SING_TEMPLATE = """
LIB "primdec.lib";
ring rA = 5,(u,v),lp;
short = 0;
ideal I = {gensA};
ideal R = std(radical(I));
"@CHARTA";
R;
ring rB = 5,(u),lp;
short = 0;
ideal J = {gensB};
ideal S = std(radical(J));
"@CHARTB";
S;
quit;
"""


def subst_w1(poly):
    out = {}
    for (a, b, c), coeff in poly.items():
        out[(a, b)] = (out.get((a, b), 0) + coeff) % P
    return {k: v for k, v in out.items() if v}


def subst_v1w0(poly):
    out = {}
    for (a, b, c), coeff in poly.items():
        if c == 0:
            out[a] = (out.get(a, 0) + coeff) % P
    return {k: v for k, v in out.items() if v}


def poly2_str(poly2):
    if not poly2:
        return "0"
    parts = []
    for (a, b), c in sorted(poly2.items(), reverse=True):
        term = str(c)
        if a:
            term += f"*u^{a}"
        if b:
            term += f"*v^{b}"
        parts.append(term)
    return "+".join(parts)


def poly1_str(poly1):
    if not poly1:
        return "0"
    return "+".join(f"{c}*u^{a}" if a else str(c)
                    for a, c in sorted(poly1.items(), reverse=True))


def parse_sing_poly(text):
    text = text.replace("-", "+-").replace(" ", "")
    out = {}
    for term in text.split("+"):
        if not term:
            continue
        coeff, eu, ev = 1, 0, 0
        if term.startswith("-"):
            coeff = -1
            term = term[1:]
        for factor in term.split("*"):
            if not factor:
                continue
            if factor.startswith("u"):
                eu = int(factor[2:]) if "^" in factor else 1
            elif factor.startswith("v"):
                ev = int(factor[2:]) if "^" in factor else 1
            else:
                coeff *= int(factor)
        key = (eu, ev)
        out[key] = (out.get(key, 0) + coeff) % P
    return {k: v for k, v in out.items() if v}


def chart_ideals(minors):
    gensA = ",\n  ".join(poly2_str(subst_w1(m)) for m in minors if m)
    polysB = [subst_v1w0(m) for m in minors]
    gensB = ",\n  ".join(poly1_str(p) for p in polysB if p) or "1"
    script = SING_TEMPLATE.format(gensA=gensA, gensB=gensB)
    out = subprocess.run(["Singular", "-q"], input=script, text=True,
                         capture_output=True, timeout=300).stdout
    sections = {"A": [], "B": []}
    current = None
    for line in out.splitlines():
        line = line.strip()
        if line == "@CHARTA":
            current = "A"
        elif line == "@CHARTB":
            current = "B"
        elif current and "=" in line and line.split("=", 1)[0].strip(
                ).startswith(("R[", "S[")):
            sections[current].append(
                parse_sing_poly(line.split("=", 1)[1]))
    return sections


def find_points(minors, fields_by_degree):
    charts = chart_ideals(minors)
    found = []
    seen_orbits = []

    def register(field, d, x):
        orbit = []
        p = tuple(x)
        while p not in orbit:
            orbit.append(p)
            p = tuple(field.frob(c) for c in p)
        if len(orbit) != d:
            return
        key = (d, frozenset(orbit))
        if key in seen_orbits:
            return
        seen_orbits.append(key)
        found.append((d, x, field))

    gensA = charts["A"]
    if gensA and not (len(gensA) == 1 and gensA[0] == {(0, 0): 1}):
        elim = next((g for g in gensA
                     if all(eu == 0 for (eu, ev) in g)), None)
        assert elim is not None, "no univariate eliminant in chart A"
        elim_coeffs = [0] * (max(ev for (_, ev) in elim) + 1)
        for (_, ev), c in elim.items():
            elim_coeffs[ev] = c
        for d in sorted(fields_by_degree):
            k = fields_by_degree[d]
            for v0 in k.elements():
                if not k.is_zero(k.poly_eval(elim_coeffs, v0)):
                    continue
                candidates = None
                for g in gensA:
                    coeffs = {}
                    for (eu, ev), c in g.items():
                        term = k.from_int(c)
                        for _ in range(ev):
                            term = k.mul(term, v0)
                        coeffs[eu] = k.add(coeffs.get(eu, k.zero), term)
                    coeffs = {e: c for e, c in coeffs.items()
                              if not k.is_zero(c)}
                    if not coeffs:
                        continue
                    roots = set()
                    for u0 in (candidates if candidates is not None
                               else k.elements()):
                        val = k.zero
                        for e, c in coeffs.items():
                            term = c
                            for _ in range(e):
                                term = k.mul(term, u0)
                            val = k.add(val, term)
                        if k.is_zero(val):
                            roots.add(u0)
                    candidates = roots
                    if not candidates:
                        break
                assert candidates is not None, \
                    "eliminant root with no u-constraint (dim > 0?)"
                for u0 in sorted(candidates):
                    register(k, d, [u0, v0, k.one])
    gensB = [{eu: c for (eu, _), c in g.items()} for g in charts["B"]]
    if gensB and not (len(gensB) == 1 and gensB[0] == {0: 1}):
        for d in sorted(fields_by_degree):
            k = fields_by_degree[d]
            for u0 in k.elements():
                ok = True
                for g in gensB:
                    val = k.zero
                    for e, c in g.items():
                        term = k.from_int(c)
                        for _ in range(e):
                            term = k.mul(term, u0)
                        val = k.add(val, term)
                    if not k.is_zero(val):
                        ok = False
                        break
                if ok:
                    register(k, d, [u0, k.one, k.zero])
    k1 = fields_by_degree[1]
    if all(m.get((4, 0, 0), 0) % P == 0 for m in minors):
        register(k1, 1, [k1.one, k1.zero, k1.zero])
    return found


# ---------------------------------------------------------------- main

def projective_points():
    pts = []
    for b in range(P):
        for c in range(P):
            pts.append((1, b, c))
    for c in range(P):
        pts.append((0, 1, c))
    pts.append((0, 0, 1))
    return pts


def fmt_elem(e):
    return list(e)


def fmt_field(field):
    return "t^%d: %s" % (field.d, field.modulus)


def main():
    rng = random.Random(SEED)
    fields_by_degree = {d: GF(d) for d in (1, 2, 3, 4, 5, 6)}
    k1 = fields_by_degree[1]
    OUTDIR.mkdir(parents=True, exist_ok=True)

    all_results = []
    for disc, path in FIELDS:
        text = path.read_text()
        mats = parse_gp_matrices(text, "secondary_norm_samples")[:6]
        basis, pairs = mats[0:3], mats[3:6]
        for lam, M in zip(SIX, mats):
            assert D_of(basis, pairs, list(lam)) == M, (disc, lam)
        Dsym = symbolic_D(basis, pairs)
        assert euler_check(Dsym), disc
        stored_T = parse_relation_matrix(text)
        if stored_T is not None:
            assert word_matrix(basis, pairs) == stored_T, disc
        minors = minors2(Dsym)

        entry = {"disc": disc,
                 "rational_points": [], "extension_points": []}

        for lam in projective_points():
            x = [k1.from_int(c) for c in lam]
            cert = certificate_at(k1, basis, pairs, x, rng=rng)
            if cert["rank"] != 1:
                continue
            td = tangent_dim(k1, minors, x)
            assert (td == 0) == cert["transverse"], (disc, lam)
            entry["rational_points"].append({
                "x": list(lam), "rank": 1,
                "kernel": [[fmt_elem(c) for c in e]
                           for e in cert["kernel"]],
                "image": [fmt_elem(c) for c in cert["image"]],
                "vs": [[fmt_elem(c) for c in v] for v in cert["vs"]],
                "y": [fmt_elem(c) for c in cert["y"]],
                "B": [[fmt_elem(c) for c in row] for row in cert["B"]],
                "detB": fmt_elem(cert["detB"]),
                "transverse": cert["transverse"],
                "tangent_dim": td,
            })

        for d, x, k in find_points(minors, fields_by_degree):
            if d == 1:
                continue
            cert = certificate_at(k, basis, pairs, x)
            rec = {
                "degree": d, "field": fmt_field(k),
                "x": [fmt_elem(c) for c in x], "rank": cert["rank"],
            }
            if cert["rank"] == 1:
                td = tangent_dim(k, minors, x)
                assert (td == 0) == cert["transverse"], (disc, d)
                rec.update({
                    "kernel": [[fmt_elem(c) for c in e]
                               for e in cert["kernel"]],
                    "image": [fmt_elem(c) for c in cert["image"]],
                    "y": [fmt_elem(c) for c in cert["y"]],
                    "B": [[fmt_elem(c) for c in row]
                          for row in cert["B"]],
                    "detB": fmt_elem(cert["detB"]),
                    "transverse": cert["transverse"],
                    "tangent_dim": td,
                })
            entry["extension_points"].append(rec)

        rat_trans = [p for p in entry["rational_points"]
                     if p["transverse"]]
        ext_trans = [p for p in entry["extension_points"]
                     if p.get("transverse")]
        entry["minimal_degree"] = (
            1 if rat_trans else
            (min(p["degree"] for p in ext_trans) if ext_trans else None))
        entry["criterion_proves_mild"] = entry["minimal_degree"] is not None
        all_results.append(entry)
        print(f"{disc}: rational rank-1 {len(entry['rational_points'])}, "
              f"transverse {len(rat_trans)}, extension transverse "
              f"{[p['degree'] for p in ext_trans]}, "
              f"mindeg {entry['minimal_degree']}")

    (OUTDIR / "certificates.json").write_text(
        json.dumps(all_results, indent=1))

    lines = [f"Transverse rank-one certificates for {len(all_results)} "
             f"computed fields", "=" * 66, "",
             "D_K          #rat rk1  #rat transverse  min k-degree  "
             "criterion applies", "-" * 66]
    for e in all_results:
        lines.append(
            f"{e['disc']:<12} {len(e['rational_points']):>8} "
            f"{len([p for p in e['rational_points'] if p['transverse']]):>15}"
            f"  {str(e['minimal_degree']):>11}  "
            f"{'YES' if e['criterion_proves_mild'] else 'NO'}")
    lines.append("")
    for e in all_results:
        lines.append(f"== {e['disc']}")
        for p in e["rational_points"]:
            lines.append(
                f"  rational x={tuple(p['x'])}  det B = {p['detB'][0]}  "
                f"{'TRANSVERSE' if p['transverse'] else 'not transverse'}"
                f"  (tangent dim {p['tangent_dim']})")
            lines.append(f"    B = {p['B']}")
        for p in e["extension_points"]:
            if p["rank"] != 1:
                lines.append(f"  degree-{p['degree']} point: rank "
                             f"{p['rank']} (criterion not applicable)")
                continue
            lines.append(
                f"  degree-{p['degree']} point over F_5[t]/({p['field']})")
            lines.append(f"    x = {p['x']}")
            lines.append(
                f"    det B = {p['detB']}  "
                f"{'TRANSVERSE' if p['transverse'] else 'not transverse'}"
                f"  (tangent dim {p['tangent_dim']})")
            lines.append(f"    B = {p['B']}")
        lines.append("")
    (OUTDIR / "report.txt").write_text("\n".join(lines) + "\n")
    print(f"written: {OUTDIR}/certificates.json, {OUTDIR}/report.txt")


if __name__ == "__main__":
    main()
