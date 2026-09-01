#!/usr/bin/env python3
"""Flag search for strong freeness of d cubic relations in d variables over
F_q, driven by Singular/Letterplace (twostd, Dp, x(1) > ... > x(d)).

Input: a JSON file {"prime": p, "d": d, "tensor": [[...d^3 entries...] x d]}
with the coefficient of X_{i1} X_{i2} X_{i3} at position d^2 (i1-1) + d (i2-1) + (i3-1).
Flags: one matrix g per coset gB of GL_d(F_q) (B upper triangular; column j
of g spans the j-th step of the flag), one representative per Frobenius
orbit when q is not prime.  Transformation: phi(X_i) = sum_a g[i][a] X_a.
Verdict rules as in fahnen_singular.py, with the target series
1/(1 - d z + d z^3) and a prefix automaton over d letters.
"""
import argparse, itertools, json, re, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import finite_field_q as fq  # noqa: E402  (FiniteField)


# ---------- words, counts, target ----------
def target_series(d, n):
    b = [1, d, d * d]
    while len(b) <= n:
        b.append(d * b[-1] - d * b[-3])
    return b[: n + 1]


def normal_word_counts(heads, d, n):
    """Counts of words over d letters with no factor in heads, degrees 0..n,
    via the prefix (Aho-Corasick) automaton; returns (counts, states)."""
    letters = "".join(str(i + 1) for i in range(d))
    prefixes = {""}
    for w in heads:
        for i in range(1, len(w)):
            prefixes.add(w[:i])
    P = sorted(prefixes); idx = {p: i for i, p in enumerate(P)}
    def step(p, ch):
        t = p + ch
        if any(t.endswith(w) for w in heads):
            return None
        while t and t not in prefixes:
            t = t[1:]
        return idx[t]
    trans = [[step(p, ch) for ch in letters] for p in P]
    cnt = [0] * len(P); cnt[idx[""]] = 1; out = [1]
    for _ in range(n):
        nxt = [0] * len(P)
        for i, c in enumerate(cnt):
            if c:
                for z in trans[i]:
                    if z is not None:
                        nxt[z] += c
        cnt = nxt; out.append(sum(cnt))
    return out, len(P)


# ---------- flags ----------
def det(M, F):
    """Determinant over F by Gaussian elimination; returns 0 iff singular."""
    n = len(M); A = [row[:] for row in M]; result = F.one
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c]), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]; result = F.neg[result]
        result = F.mul(result, A[c][c]); inv = F.inv[A[c][c]]
        for r in range(c + 1, n):
            if A[r][c]:
                f = F.mul(A[r][c], inv)
                A[r] = [F.sub(x, F.mul(f, y)) for x, y in zip(A[r], A[c])]
    return result


def rref(rows, F, d):
    m = [list(r) for r in rows]; r = 0
    for c in range(d):
        piv = next((i for i in range(r, len(m)) if m[i][c]), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        inv = F.inv[m[r][c]]; m[r] = [F.mul(x, inv) for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c]:
                f = m[i][c]; m[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(m[i], m[r])]
        r += 1
    return tuple(tuple(row) for row in m)


def in_span(v, cols, F, d):
    if not cols:
        return not any(v)
    return rref(list(cols) + [v], F, d)[len(cols)] == tuple([0] * d) if len(cols) < d else True


def flag_representatives(F, d, frobenius_orbits=True):
    """Complete flags as tuples of RREF subspace ids; representative matrix per flag."""
    q = F.q
    vectors = [v for v in itertools.product(range(q), repeat=d) if any(v)]
    flags = {}  # flag id (tuple of rref of the first j columns, j=1..d-1) -> g (columns)
    def extend(cols, ids):
        j = len(cols)
        if j == d - 1:
            for u in vectors:
                g = cols + [u]
                if det([[g[c][r] for c in range(d)] for r in range(d)], F):
                    flags[tuple(ids)] = g
                    return
            raise RuntimeError("no completion")
        seen = set()
        for w in vectors:
            if in_span(w, cols, F, d):
                continue
            sid = rref(cols + [w], F, d)
            if sid in seen:
                continue
            seen.add(sid)
            extend(cols + [w], ids + [sid])
    extend([], [])
    expected = 1
    for k in range(1, d + 1):
        expected *= (q ** k - 1) // (q - 1)
    assert len(flags) == expected, (len(flags), expected)
    reps = []
    if not frobenius_orbits or F.e == 1:
        for fid in sorted(flags):
            reps.append(flags[fid])
        return reps
    done = set()
    for fid in sorted(flags):
        if fid in done:
            continue
        cur = fid
        for _ in range(F.e):
            done.add(cur)
            cur = tuple(rref([tuple(F.frobenius(x) for x in row) for row in sid], F, d) for sid in cur)
        reps.append(flags[fid])
    return reps


# ---------- transformation and Singular ----------
def transform(T, g, F, d):
    """g given as list of columns (each a d-vector): g[i][a] = column a, entry i."""
    G = [[g[a][i] for a in range(d)] for i in range(d)]  # G[i][a]
    out = []
    for row in T:
        new = [0] * (d ** 3)
        for idx, c in enumerate(row):
            if not c:
                continue
            i, j, k = idx // (d * d), (idx // d) % d, idx % d
            for a in range(d):
                ga = G[i][a]
                if not ga:
                    continue
                cga = F.mul(c, ga)
                for b in range(d):
                    gb = G[j][b]
                    if not gb:
                        continue
                    cgab = F.mul(cga, gb)
                    for e in range(d):
                        ge = G[k][e]
                        if ge:
                            pos = d * d * a + d * b + e
                            new[pos] = F.add(new[pos], F.mul(cgab, ge))
        out.append(new)
    return out


def coeff_text(x, F):
    if F.e == 1:
        return str(x)
    terms = [str(c) if i == 0 else (f"{c}*a" if i == 1 else f"{c}*a^{i}") for i, c in enumerate(F.digits[x]) if c]
    return "(" + "+".join(terms) + ")"


def minpoly_text(F):
    terms = [f"a^{F.e}"] + [str(c) if i == 0 else (f"{c}*a" if i == 1 else f"{c}*a^{i}") for i, c in enumerate(F.modulus[:-1]) if c]
    return "+".join(terms)


def relation_text(row, F, d):
    terms = []
    for idx, c in enumerate(row):
        if c:
            i, j, k = idx // (d * d), (idx // d) % d, idx % d
            terms.append(f"{coeff_text(c, F)}*x({i+1})*x({j+1})*x({k+1})")
    return "+".join(terms) or "0"


def script(batch, T, F, d, bound):
    ring = (f"ring r=({F.p},a),(x(1..{d})),Dp;\nminpoly={minpoly_text(F)};" if F.e > 1 else f"ring r={F.p},(x(1..{d})),Dp;")
    lines = ['LIB "freegb.lib";', ring, f"def R=freeAlgebra(r,{bound});", "setring R;", "int q;"]
    for n, g in batch:
        rels = [relation_text(row, F, d) for row in transform(T, g, F, d)]
        lines += ["ideal I=" + ",\n ".join(rels) + ";", "ideal G=twostd(I);", f'print("FLAG {n} BEGIN");',
                  "for (q=1; q<=size(G); q++) { lead(G[q]); }", 'print("FLAG END");', "kill I; kill G;"]
    lines.append("quit;")
    return "\n".join(lines) + "\n"


def parse(out):
    result = {}
    for m in re.finditer(r"FLAG (\d+) BEGIN\n(.*?)FLAG END", out, re.S):
        leads = ["".join(re.findall(r"x\((\d+)\)", l)) for l in m.group(2).splitlines()]
        leads = [w for w in leads if w]
        minimal = []
        for w in sorted(set(leads), key=lambda w: (len(w), w)):
            if not any(v in w for v in minimal):
                minimal.append(w)
        result[int(m.group(1))] = minimal
    return result


def verdict(heads, d, bound):
    m = max(len(w) for w in heads)
    terminated = 2 * m - 1 <= bound
    series, states = normal_word_counts(heads, d, bound)
    target = target_series(d, bound)
    dev = next((n for n in range(bound + 1) if series[n] != target[n]), None)
    if dev is not None:
        return {"verdict": "NOT_STRONGLY_FREE", "deviation_degree": dev, "terminated": terminated, "states": states}
    if terminated:
        w = 2 * (states + 3) + 3
        s2, _ = normal_word_counts(heads, d, w); t2 = target_series(d, w)
        ok = s2 == t2
        return {"verdict": "STRONGLY_FREE" if ok else "NOT_STRONGLY_FREE",
                "deviation_degree": None if ok else next(i for i, (x, y) in enumerate(zip(s2, t2)) if x != y),
                "terminated": True, "states": states}
    return {"verdict": "INCONCLUSIVE", "deviation_degree": None, "terminated": False, "states": states}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tensor_json")
    ap.add_argument("--degree", type=int, default=1, help="extension degree e of F_q over F_p")
    ap.add_argument("--bound", type=int, default=9)
    ap.add_argument("--batch", type=int, default=10)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--indices", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--timeout", type=int, default=7200)
    ap.add_argument("--singular", default="Singular")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    src = json.loads(Path(a.tensor_json).read_text())
    p, d = src["prime"], src["d"]
    F = fq.FiniteField(p, a.degree)
    T = [[int(v) % p for v in row] for row in src["tensor"]]
    assert len(T) == d and all(len(r) == d ** 3 for r in T)
    t0 = time.perf_counter()
    reps = flag_representatives(F, d)
    chosen = [int(i) for i in a.indices.split(",") if i.strip()] if a.indices else list(range(len(reps)))
    chosen = [i for i in chosen if i % a.shards == a.shard]
    if a.limit:
        chosen = chosen[: a.limit]
    print(f"d={d}, q={F.q}: {len(reps)} flags ({time.perf_counter()-t0:.1f}s), {len(chosen)} in this shard, bound {a.bound}", flush=True)
    ledger, stop = [], False
    a.out.parent.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(chosen), a.batch):
        batch = [(n, reps[n]) for n in chosen[start:start + a.batch]]
        t0 = time.perf_counter()
        run = subprocess.run([a.singular, "-q"], input=script(batch, T, F, d, a.bound), text=True,
                             capture_output=True, timeout=a.timeout)
        secs = time.perf_counter() - t0
        found = parse(run.stdout)
        for n, g in batch:
            G = [[g[c][r] for c in range(d)] for r in range(d)]  # row convention matrix
            if n not in found:
                entry = {"index": n, "basis": G, "verdict": "SINGULAR_ERROR", "detail": (run.stderr + run.stdout)[-400:]}
            else:
                heads = found[n]; v = verdict(heads, d, a.bound)
                entry = {"index": n, "basis": G, "head_words": heads, "max_head_degree": max(map(len, heads)), "bound": a.bound, **v}
            ledger.append(entry)
            print(f"[{n}/{len(reps)}] {entry['verdict']} |heads|={len(entry.get('head_words', []))} m={entry.get('max_head_degree')} states={entry.get('states')}", flush=True)
            if entry["verdict"] in ("STRONGLY_FREE", "NOT_STRONGLY_FREE"):
                print(f"{entry['verdict']} — basis {G} head words {entry.get('head_words')}", flush=True); stop = True
        print(f"  batch of {len(batch)} in {secs:.1f}s ({secs/len(batch):.2f}s per flag)", flush=True)
        a.out.write_text(json.dumps({"prime": p, "d": d, "degree": F.e, "q": F.q, "modulus_low_to_high": F.modulus, "bound": a.bound,
                                     "engine": "Singular freegb/Letterplace twostd, Dp", "source": str(a.tensor_json),
                                     "shards": a.shards, "shard": a.shard, "ledger": ledger}, indent=1))
        if stop:
            break


if __name__ == "__main__":
    main()
