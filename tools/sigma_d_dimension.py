#!/usr/bin/env python3
"""Krull dimension of the norm-degeneracy scheme for every computed field.

The paper states that Sigma_D is zero-dimensional for every field of
the census.  This script proves the statement computationally, one
field at a time: it builds the matrix family D_x symbolically over
F_p, forms the nine 2x2 minors (the equations of Sigma_D), computes a
Groebner basis of the minor ideal in F_p[x1,x2,x3], and reads off the
Krull dimension of the affine cone from the leading terms.  The claim
is equivalent to dimension at most one for every field: dimension 0
means Sigma_D is empty, dimension 1 means finitely many closed
points, and dimension 2 would exhibit a curve inside Sigma_D.

Two routes construct the family:

  * tensor route (p = 3): read tensor.json from
    records/p3/source-tensors/ and contract, exactly as in
    tools/p3_cone_criterion.py;
  * matrix route (any p): obtain the six sampled matrices of a
    field --- from matrices.tsv of a certificate directory (--tsv),
    or for p = 5 by the certificate export of
    tools/transverse_rank_one.py (--p5) --- and solve the
    polarization system for the six coefficient matrices of the
    quadratic family; when the stored characters overdetermine the
    system, the script asserts the consistency of every stored
    matrix with the solved family.

The two routes are cross-validated against each other on p = 3
fields (--validate), where both data files exist.

Usage:
  python3 tools/sigma_d_dimension.py --p3 [--workers N] [--out FILE]
  python3 tools/sigma_d_dimension.py --p5 [--workers N] [--out FILE]
  python3 tools/sigma_d_dimension.py --tsv certificates/p7/K-931506071-p7 --prime 7
  python3 tools/sigma_d_dimension.py --validate 6
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "records" / "p3" / "source-tensors"

P = 3
PAIR_BUDGET = 200000

# ----------------------------------------------------------------------
# Polynomials as {(a, b, c): coefficient mod P}; graded reverse
# lexicographic order.


def grevlex_greater(m, n):
    dm, dn = sum(m), sum(n)
    if dm != dn:
        return dm > dn
    for a, b in zip(reversed(m), reversed(n)):
        if a != b:
            return a < b
    return False


def leading_monomial(f):
    lead = None
    for m in f:
        if lead is None or grevlex_greater(m, lead):
            lead = m
    return lead


def poly_add(f, g):
    h = dict(f)
    for m, c in g.items():
        c2 = (h.get(m, 0) + c) % P
        if c2:
            h[m] = c2
        else:
            h.pop(m, None)
    return h


def poly_scale(f, c):
    c %= P
    return {m: (c * v) % P for m, v in f.items()} if c else {}


def mono_mul(m, n):
    return (m[0] + n[0], m[1] + n[1], m[2] + n[2])


def poly_mul_mono(f, m, c):
    c %= P
    return {mono_mul(mm, m): (c * v) % P for mm, v in f.items()} if c else {}


def poly_mul(f, g):
    h = {}
    for m, c in f.items():
        for n, d in g.items():
            mn = mono_mul(m, n)
            h[mn] = (h.get(mn, 0) + c * d) % P
    return {m: c for m, c in h.items() if c}


def mono_divides(m, n):
    return all(a <= b for a, b in zip(m, n))


def mono_div(n, m):
    return (n[0] - m[0], n[1] - m[1], n[2] - m[2])


def inv_mod_p(c):
    return pow(c, P - 2, P)


def reduce_poly(f, basis):
    f = dict(f)
    result = {}
    while f:
        lm = leading_monomial(f)
        lc = f[lm]
        for lt, ltc, g in basis:
            if mono_divides(lt, lm):
                factor = mono_div(lm, lt)
                coeff = (-lc * inv_mod_p(ltc)) % P
                f = poly_add(f, poly_mul_mono(g, factor, coeff))
                break
        else:
            result[lm] = lc
            del f[lm]
    return result


def buchberger(gens):
    basis = []
    for g in gens:
        g = {m: c % P for m, c in g.items() if c % P}
        if g:
            basis.append((leading_monomial(g), g[leading_monomial(g)], g))
    pairs = [(i, j) for i in range(len(basis)) for j in range(i)]
    steps = 0
    while pairs:
        steps += 1
        if steps > PAIR_BUDGET:
            return None
        i, j = pairs.pop()
        lti, ltci, gi = basis[i]
        ltj, ltcj, gj = basis[j]
        lcm = tuple(max(a, b) for a, b in zip(lti, ltj))
        if lcm == mono_mul(lti, ltj):
            continue
        s = poly_add(
            poly_mul_mono(gi, mono_div(lcm, lti), inv_mod_p(ltci)),
            poly_mul_mono(gj, mono_div(lcm, ltj), (-inv_mod_p(ltcj)) % P),
        )
        s = reduce_poly(s, basis)
        if s:
            lt = leading_monomial(s)
            basis.append((lt, s[lt], s))
            pairs.extend((len(basis) - 1, k) for k in range(len(basis) - 1))
    return basis


def dimension_from_basis(basis):
    """dim V(I) inside A^3, from the leading monomials."""
    leads = [lt for lt, _, _ in basis]
    if any(sum(lt) == 0 for lt in leads):
        return -1
    best = 0
    for mask in range(1, 8):
        s = {v for v in range(3) if mask & (1 << v)}
        if any(all(lt[v] == 0 for v in range(3) if v not in s) for lt in leads):
            continue
        best = max(best, len(s))
    return best


# ----------------------------------------------------------------------
# Route 1: the family from a stored tensor (p = 3), as in
# tools/p3_cone_criterion.py.


def entries_from_tensor(T):
    unit = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    entries = [[{} for _ in range(3)] for _ in range(3)]
    for j in range(3):
        for ell in range(3):
            f = {}
            for i in range(3):
                for k in range(3):
                    c = T[ell][9 * i + 3 * j + k] % P
                    if c:
                        m = mono_mul(unit[i], unit[k])
                        f[m] = (f.get(m, 0) + c) % P
            entries[j][ell] = {m: c for m, c in f.items() if c}
    return entries


# ----------------------------------------------------------------------
# Route 2: the family from matrices.tsv by solving the polarization
# system.  The quadratic family is sum_{i<=k} x_i x_k C_{ik} with six
# unknown 3x3 matrices; every stored character contributes one linear
# equation per matrix entry, and the system is overdetermined.

MONOMIALS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def entries_from_tsv(path):
    rows = [line.split("\t") for line in path.read_text().splitlines()[1:] if line]
    samples = {}
    for row in rows:
        x = (int(row[1]) % P, int(row[2]) % P, int(row[3]) % P)
        ell = int(row[5]) - 1
        d = (int(row[6]) % P, int(row[7]) % P, int(row[8]) % P)
        samples.setdefault(x, {})[ell] = d
    return entries_from_samples(samples)


def entries_from_samples(samples):
    for x, cols in samples.items():
        if sorted(cols) != [0, 1, 2]:
            raise ValueError(f"unvollstaendige Spalten fuer {x}")
    xs = sorted(samples)
    matrix_rows = [[(x[i] * x[k]) % P for (i, k) in MONOMIALS] for x in xs]
    entries = [[None] * 3 for _ in range(3)]
    for j in range(3):
        for ell in range(3):
            rhs = [samples[x][ell][j] for x in xs]
            solution = solve_mod_p(matrix_rows, rhs)
            f = {}
            for (i, k), c in zip(MONOMIALS, solution):
                if c:
                    m = [0, 0, 0]
                    m[i] += 1
                    m[k] += 1
                    f[tuple(m)] = c
            entries[j][ell] = f
    return entries


def solve_mod_p(rows, rhs):
    """Solve the overdetermined system exactly; raise on inconsistency."""
    aug = [[v % P for v in row] + [b % P] for row, b in zip(rows, rhs)]
    n_cols = len(rows[0])
    pivots = []
    r = 0
    for c in range(n_cols):
        pivot = next((k for k in range(r, len(aug)) if aug[k][c]), None)
        if pivot is None:
            continue
        aug[r], aug[pivot] = aug[pivot], aug[r]
        inv = inv_mod_p(aug[r][c])
        aug[r] = [(inv * v) % P for v in aug[r]]
        for k in range(len(aug)):
            if k != r and aug[k][c]:
                factor = aug[k][c]
                aug[k] = [(v - factor * w) % P for v, w in zip(aug[k], aug[r])]
        pivots.append(c)
        r += 1
    if r < n_cols:
        raise ValueError("Polarisationssystem unterbestimmt")
    for k in range(r, len(aug)):
        if aug[k][n_cols]:
            raise ValueError("gespeicherte Matrizen inkonsistent mit der Familie")
    solution = [0] * n_cols
    for row_index, c in enumerate(pivots):
        solution[c] = aug[row_index][n_cols]
    return solution


# ----------------------------------------------------------------------
# The check itself.


def poly_to_singular(f):
    terms = []
    for (a, b, c), coeff in sorted(f.items()):
        factors = [str(coeff)]
        for name, exponent in (("x1", a), ("x2", b), ("x3", c)):
            if exponent:
                factors.append(f"{name}^{exponent}")
        terms.append("*".join(factors))
    return "+".join(terms) if terms else "0"


def singular_dimensions(named_entries):
    """Recompute every dimension with Singular; one process, many rings."""
    import subprocess
    blocks = []
    for name, entries in named_entries:
        gens = ",".join(poly_to_singular(f) for f in minors(entries))
        blocks.append(
            f'ring r = {P},(x1,x2,x3),dp; ideal I = {gens};\n'
            f'print("{name}:"+string(dim(groebner(I)))); kill r;')
    script = "\n".join(blocks) + "\nquit;\n"
    run = subprocess.run(["Singular", "-q"], input=script,
                         capture_output=True, text=True, timeout=7200)
    if run.returncode != 0:
        raise SystemExit(f"Singular failed:\n{run.stderr[-500:]}")
    out = {}
    for line in run.stdout.splitlines():
        name, _, value = line.partition(":")
        if value.lstrip("-").isdigit():
            out[name] = int(value)
    return out


def minors(entries):
    out = []
    for r1 in range(3):
        for r2 in range(r1 + 1, 3):
            for c1 in range(3):
                for c2 in range(c1 + 1, 3):
                    det = poly_add(
                        poly_mul(entries[r1][c1], entries[r2][c2]),
                        poly_scale(poly_mul(entries[r1][c2], entries[r2][c1]), P - 1),
                    )
                    if det:
                        out.append(det)
    return out


def check_entries(name, entries):
    started = time.monotonic()
    basis = buchberger(minors(entries))
    if basis is None:
        return {"field": name, "dim": None, "status": "PAIR_BUDGET"}
    dim = dimension_from_basis(basis)
    return {"field": name, "dim": dim,
            "status": "OK" if dim <= 1 else "DIM_GE_2",
            "ms": round(1000 * (time.monotonic() - started))}


def check_tensor_dir(path):
    T = json.loads((path / "tensor.json").read_text())["tensor_3_by_27"]
    return check_entries(path.name, entries_from_tensor(T))


def check_p5_certificate(job):
    disc, cert_path = job
    import transverse_rank_one as tr
    mats = tr.certificate_matrices(cert_path)
    samples = {tuple(x): {ell: tuple(mats[q][j][ell] for j in range(3))
                          for ell in range(3)}
               for q, x in enumerate(tr.SIX)}
    return check_entries(f"D{disc}", entries_from_samples(samples))


def main():
    global P
    ap = argparse.ArgumentParser()
    ap.add_argument("--p3", action="store_true", help="sweep records/p3/source-tensors")
    ap.add_argument("--p5", action="store_true",
                    help="sweep the p5 certificates via the GP matrix export")
    ap.add_argument("--tsv", type=Path, help="certificate directory with matrices.tsv")
    ap.add_argument("--prime", type=int, default=3)
    ap.add_argument("--validate", type=int, default=0,
                    help="cross-validate both routes on N p=3 fields")
    ap.add_argument("--crosscheck", choices=["p3", "p5", "p7"],
                    help="recompute every dimension of a collection with "
                         "Singular and compare with the banked record")
    ap.add_argument("--banked", type=Path,
                    help="banked sigma-dimension.json for --crosscheck")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    P = args.prime

    if args.validate:
        P = 3
        dirs = [d for d in sorted(SRC.glob("*/D-*"))
                if (d / "matrices.tsv").exists()][: args.validate]
        for d in dirs:
            T = json.loads((d / "tensor.json").read_text())["tensor_3_by_27"]
            a = entries_from_tensor(T)
            b = entries_from_tsv(d / "matrices.tsv")
            same = a == b
            print(f"{d.name}: Tensor-Route == TSV-Route: {same}")
            if not same:
                return 1
        print("Beide Routen identisch.")
        return 0

    if args.crosscheck:
        import sys
        named = []
        if args.crosscheck == "p3":
            P = 3
            for d in sorted(SRC.glob("*/D-*")):
                T = json.loads((d / "tensor.json").read_text())["tensor_3_by_27"]
                named.append((d.name, entries_from_tensor(T)))
        elif args.crosscheck == "p5":
            sys.path.insert(0, str(HERE))
            import transverse_rank_one as tr
            P = 5
            for disc, _, cert in tr.discover_fields():
                mats = tr.certificate_matrices(cert)
                samples = {tuple(x): {ell: tuple(mats[q][j][ell] for j in range(3))
                                      for ell in range(3)}
                           for q, x in enumerate(tr.SIX)}
                named.append((f"D{disc}", entries_from_samples(samples)))
        else:
            P = 7
            for d in sorted((ROOT / "certificates" / "p7").glob("K-*-p7")):
                named.append((d.name, entries_from_tsv(d / "matrices.tsv")))
        print(f"{len(named)} Koerper, eine Singular-Sitzung", flush=True)
        singular = singular_dimensions(named)
        banked = {row["field"]: row["dim"]
                  for row in json.loads(args.banked.read_text())["fields"]}
        disagreements = sorted(
            name for name in banked
            if singular.get(name) != banked[name])
        missing = sorted(set(banked) - set(singular))
        payload = {"engine": "Singular, dim(groebner(I)) over F_p",
                   "total": len(banked),
                   "agree": len(banked) - len(disagreements),
                   "disagreements": disagreements,
                   "missing_from_singular_output": missing}
        if args.out:
            args.out.write_text(json.dumps(payload, indent=1) + "\n")
        print(json.dumps(payload))
        return 0 if not disagreements and not missing else 1

    if args.tsv:
        entries = entries_from_tsv(args.tsv / "matrices.tsv")
        row = check_entries(args.tsv.name, entries)
        print(json.dumps(row))
        if args.out:
            args.out.write_text(json.dumps(row, indent=1) + "\n")
        return 0 if row["status"] == "OK" else 1

    if args.p5:
        import sys
        sys.path.insert(0, str(HERE))
        import transverse_rank_one as tr
        P = 5
        jobs = [(disc, cert) for disc, _, cert in tr.discover_fields()]
        print(f"{len(jobs)} Koerper, {args.workers} Worker", flush=True)
        results = []
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(check_p5_certificate, j) for j in jobs]
            for done, fut in enumerate(as_completed(futures), start=1):
                row = fut.result()
                results.append(row)
                if row["status"] != "OK":
                    print(f"AUFFAELLIG: {row}", flush=True)
                if done % 50 == 0 or done == len(jobs):
                    print(f"{done}/{len(jobs)}", flush=True)
        dist = {}
        for row in results:
            dist[str(row["dim"])] = dist.get(str(row["dim"]), 0) + 1
        bad = [r for r in results if r["status"] != "OK"]
        payload = {"summary": {"total": len(results), "dim_distribution": dist,
                               "non_ok": bad, "claim_holds": not bad},
                   "fields": sorted(results, key=lambda r: r["field"])}
        if args.out:
            args.out.write_text(json.dumps(payload, indent=1))
        print(json.dumps(payload["summary"]))
        return 0 if not bad else 1

    if args.p3:
        dirs = sorted(SRC.glob("*/D-*"))
        print(f"{len(dirs)} Koerper, {args.workers} Worker", flush=True)
        results = []
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(check_tensor_dir, d) for d in dirs]
            for done, fut in enumerate(as_completed(futures), start=1):
                row = fut.result()
                results.append(row)
                if row["status"] != "OK":
                    print(f"AUFFAELLIG: {row}", flush=True)
                if done % 500 == 0 or done == len(dirs):
                    print(f"{done}/{len(dirs)}", flush=True)
        dist = {}
        for row in results:
            dist[str(row["dim"])] = dist.get(str(row["dim"]), 0) + 1
        bad = [r for r in results if r["status"] != "OK"]
        payload = {"summary": {"total": len(results), "dim_distribution": dist,
                               "non_ok": bad, "claim_holds": not bad},
                   "fields": sorted(results, key=lambda r: r["field"])}
        if args.out:
            args.out.write_text(json.dumps(payload, indent=1))
        print(json.dumps(payload["summary"]))
        return 0 if not bad else 1

    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
