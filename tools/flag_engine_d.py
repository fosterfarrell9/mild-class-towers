#!/usr/bin/env python3
"""d-generic completion engine over F_q (second engine, plain Python):
words over the letters 1..d, order X1 > X2 > ... > Xd refining degree,
overlap-based completion as in fahnen_fq.py / the paper's Algorithm A.6."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import finite_field_q as fq
import flag_search_d as fd


def key(word):
    return len(word), [-int(c) for c in word]


def head(poly):
    return max(poly, key=key)


def sub_scaled(target, factor, source, F):
    for w, c in source.items():
        v = F.sub(target.get(w, 0), F.mul(factor, c))
        if v:
            target[w] = v
        else:
            target.pop(w, None)


def reduce_poly(poly, basis, F):
    poly = dict(poly)
    changed = True
    while changed and poly:
        changed = False
        for w in sorted(poly, key=key, reverse=True):
            c = poly.get(w)
            if c is None:
                continue
            for lead, g in basis.items():
                pos = w.find(lead)
                if pos >= 0:
                    pre, suf = w[:pos], w[pos + len(lead):]
                    sub_scaled(poly, c, {pre + u + suf: cu for u, cu in g.items()}, F)
                    changed = True
                    break
            if changed:
                break
    return poly


def monic(poly, F):
    inv = F.inv[poly[head(poly)]]
    return {w: F.mul(c, inv) for w, c in poly.items()}


def overlaps(left, right):
    for length in range(1, min(len(left), len(right))):
        if left[-length:] == right[:length]:
            yield left + right[length:]


def complete(relations, F, max_degree):
    basis = {}
    pending = [dict(r) for r in relations if r]

    def insert_all():
        nonlocal pending
        changed = True
        while changed:
            changed = False
            queue, pending = pending, []
            for poly in queue:
                red = reduce_poly(poly, basis, F)
                if not red:
                    continue
                red = monic(red, F)
                lead = head(red)
                retired = [old for old in basis if old != lead and lead in old]
                basis[lead] = red
                for old in retired:
                    pending.append(basis.pop(old))
                changed = True

    insert_all()
    degree = 3
    while degree < max_degree:
        degree += 1
        new = []
        leads = sorted(basis, key=key)
        for left in leads:
            for right in leads:
                for comp in overlaps(left, right):
                    if len(comp) != degree:
                        continue
                    g1, g2 = basis.get(left), basis.get(right)
                    if g1 is None or g2 is None:
                        continue
                    suffix = comp[len(left):]; prefix = comp[:len(comp) - len(right)]
                    s = {w + suffix: c for w, c in g1.items()}
                    sub_scaled(s, 1, {prefix + w: c for w, c in g2.items()}, F)
                    rem = reduce_poly(s, basis, F)
                    if rem:
                        new.append(rem)
        if new:
            pending.extend(new)
            insert_all()
        largest = max(len(w) for w in basis)
        if degree >= 2 * largest - 1:
            return basis, degree, True
    return basis, max_degree, False


def verdict(T, G, F, d, max_degree):
    """G: row-convention matrix (G[i][a]); T: d rows of d^3 entries."""
    g_cols = [[G[i][a] for i in range(d)] for a in range(d)]  # columns for fd.transform
    rows = fd.transform(T, g_cols, F, d)
    rels = []
    for row in rows:
        poly = {}
        for idx, c in enumerate(row):
            if c:
                i, j, k = idx // (d * d), (idx // d) % d, idx % d
                poly[f"{i+1}{j+1}{k+1}"] = c
        rels.append(poly)
    basis, processed, terminated = complete(rels, F, max_degree)
    leads = sorted(basis, key=lambda w: (len(w), w))
    series, states = fd.normal_word_counts(leads, d, processed)
    target = fd.target_series(d, processed)
    dev = next((n for n in range(processed + 1) if series[n] != target[n]), None)
    if dev is not None:
        return {"verdict": "NOT_STRONGLY_FREE", "deviation_degree": dev, "processed_degree": processed, "terminated": terminated, "head_words": leads}
    if terminated:
        w = 2 * (states + 3) + 3
        s2, _ = fd.normal_word_counts(leads, d, w); t2 = fd.target_series(d, w)
        ok = s2 == t2
        return {"verdict": "STRONGLY_FREE" if ok else "NOT_STRONGLY_FREE", "processed_degree": processed, "terminated": True, "head_words": leads, "states": states}
    return {"verdict": "INCONCLUSIVE", "processed_degree": processed, "terminated": False, "head_words": leads}
