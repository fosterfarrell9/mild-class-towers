#!/usr/bin/env python3
"""Standalone verifier for flag certificates of rank d, 3 <= d <= 9 (plain Python).

For a certificate in records/p3/flag-certificates-rank-four/ it checks, from
scratch, that after the change of variables over the stated field the
completion (Lemma A.5) of the d cubic relations terminates with exactly the
listed head words, and that the words avoiding them are counted by
1/(1 - d z + d z^3) in every degree (Lemma A.4, Corollary A.2).  The
relation tensor is read from records/p3/rank-four-<|D|>/tensor.json, the
tensor reconstructed by tools/verify_certificate_general.gp from the
arithmetic certificate of that directory.

  python3 tools/verify_flag_certificate_d.py records/p3/flag-certificates-rank-four/653329427.json
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import finite_field_q as fq  # noqa: E402
import flag_engine_d as engine  # noqa: E402

REPO = HERE.parent
VERIFIER_VERSION = "2026-09-01.d3"
PFLICHT = {"discriminant", "rank", "field", "basis", "head_words", "proper_prefixes", "search", "claim", "alternatives"}
PRIMZAHLEN = (3, 5)


def schema(cert, dateiname=None):
    if set(cert) != PFLICHT:
        raise ValueError(f"keys: expected {sorted(PFLICHT)}, found {sorted(cert)}")
    D = cert["discriminant"]; d = cert["rank"]
    if not isinstance(D, int) or D >= 0: raise ValueError("discriminant must be a negative integer")
    if dateiname is not None and Path(dateiname).name != f"{abs(D)}.json": raise ValueError("file name does not match the discriminant")
    if not isinstance(d, int) or not 3 <= d <= 9: raise ValueError("rank must be an integer between 3 and 9")
    fld = cert["field"]
    if set(fld) != {"q", "degree", "modulus_low_to_high"}: raise ValueError("field keys")
    e, q = fld["degree"], fld["q"]
    prime = next((pr for pr in PRIMZAHLEN if pr ** e == q), None)
    if prime is None: raise ValueError("q must equal p^degree for p in {3, 5}")
    if (prime, e) not in fq.MODULI: raise ValueError(f"no documented modulus for p={prime}, degree {e}")
    if list(fld["modulus_low_to_high"]) != list(fq.MODULI[(prime, e)]): raise ValueError("modulus differs from the documented one")
    B = cert["basis"]
    if not (isinstance(B, list) and len(B) == d and all(isinstance(r, list) and len(r) == d for r in B)):
        raise ValueError(f"basis must be a {d}x{d} matrix")
    F = fq.FiniteField(prime, e)
    G = []
    for row in B:
        grow = []
        for x in row:
            if e == 1:
                if not (isinstance(x, int) and 0 <= x < prime): raise ValueError(f"basis entries over F_{prime} must be 0..{prime-1}")
                grow.append(x)
            else:
                if not (isinstance(x, list) and len(x) == e and all(isinstance(c, int) and 0 <= c < prime for c in x)):
                    raise ValueError(f"basis entries must be coefficient lists of length {e}")
                grow.append(F.enc[tuple(x)])
        G.append(grow)
    letters = set(str(i + 1) for i in range(d))
    W = cert["head_words"]
    if not (isinstance(W, list) and W and all(isinstance(w, str) and w and set(w) <= letters for w in W)):
        raise ValueError(f"head_words must be nonempty words over 1..{d}")
    if len(set(W)) != len(W): raise ValueError("duplicate head words")
    if any(u != w and u in w for u in W for w in W): raise ValueError("head words must be factor-free")
    s = len({""} | {w[:i] for w in W for i in range(1, len(w))})
    if cert["proper_prefixes"] != s: raise ValueError(f"proper_prefixes must be {s}")
    if not isinstance(cert["claim"], str): raise ValueError("claim must be a string (it is documentation; the verifier does not interpret it)")
    return prime, e, d, F, G


def lade_tensor(D, d, prime):
    import re as _re
    p = REPO / "records" / "p3" / f"rank-four-{abs(D)}" / "tensor.json"
    if not p.exists():
        p = REPO / "records" / "p3" / "rank-four-census" / f"K-{abs(D)}" / "tensor.json"
    if not p.exists(): raise KeyError(f"tensor of {D}: {p} not found")
    data = json.loads(p.read_text())
    if data.get("discriminant") != D: raise KeyError(f"{p}: discriminant of the tensor record differs from the certificate")
    if data.get("prime") != prime: raise KeyError(f"{p}: prime of the tensor record ({data.get('prime')}) differs from the certificate ({prime})")
    if data.get("rank") != d: raise KeyError(f"{p}: rank of the tensor record differs from the certificate")
    T = [[int(v) % prime for v in row] for row in data["tensor"]]
    if len(T) != d or any(len(r) != d ** 3 for r in T): raise KeyError(f"{p}: tensor is not {d} x {d**3}")
    # Bindung an den arithmetischen Record im selben Verzeichnis:
    log = p.parent / "verification.log"
    if not log.exists(): raise KeyError(f"{p.parent}: verification.log not found")
    logtext = log.read_text()
    if "CERTIFICATE VERIFIED   (812 checks" not in logtext:
        raise KeyError(f"{log}: the witness certificate of this directory is not verified")
    m = _re.search(r"^TENSOR = (\[.*\])$", logtext, _re.M)
    if not m: raise KeyError(f"{log}: no reconstructed TENSOR")
    import ast as _ast
    V = [[int(v) % prime for v in row] for row in _ast.literal_eval(m.group(1))]
    if V != T: raise KeyError(f"{log}: reconstructed tensor differs from tensor.json")
    certgp = p.parent / "certificate.gp"
    if not certgp.exists(): raise KeyError(f"{p.parent}: certificate.gp not found")
    kopf = certgp.read_text()[:200]
    mk = _re.match(r"\s*\[\s*\d+\s*,\s*\d+\s*,\s*(\d+)\s*,\s*[^,]+,\s*(-\d+)\s*,", kopf)
    if not mk: raise KeyError(f"{certgp}: cannot read prime and discriminant from the certificate head")
    if int(mk.group(1)) != prime or int(mk.group(2)) != D:
        raise KeyError(f"{certgp}: prime/discriminant of the witness certificate differ from the flag certificate")
    return T, str(p)


def pruefe(cert, dateiname=None):
    D = cert["discriminant"]
    try:
        prime, e, d, F, G = schema(cert, dateiname)
    except ValueError as exc:
        return {"discriminant": D, "verified": False, "stage": "schema", "detail": str(exc)}
    g_cols = [[G[i][a] for i in range(d)] for a in range(d)]
    if not engine.fd.det([[g_cols[c][r] for c in range(d)] for r in range(d)], F):
        return {"discriminant": D, "rank": d, "field": prime ** e, "verified": False, "stage": "basis-singular"}
    T, quelle = lade_tensor(D, d, prime)
    t0 = time.monotonic()
    m = max(map(len, cert["head_words"]))
    res = engine.verdict(T, G, F, d, 2 * m - 1)
    ok = res["verdict"] == "STRONGLY_FREE" and sorted(res["head_words"]) == sorted(cert["head_words"])
    out = {"discriminant": D, "rank": d, "prime": prime, "field": prime ** e, "verified": ok, "verifier_version": VERIFIER_VERSION,
           "head_words": len(cert["head_words"]), "max_head_degree": m, "seconds": round(time.monotonic() - t0, 1), "tensor_source": quelle}
    if not ok:
        out["stage"] = "completion"; out["detail"] = {k: v for k, v in res.items()}
    h = hashlib.sha256(json.dumps(cert, sort_keys=True).encode()); h.update(json.dumps(T).encode()); h.update(VERIFIER_VERSION.encode())
    out["digest"] = h.hexdigest()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zertifikat")
    a = ap.parse_args()
    try:
        r = pruefe(json.load(open(a.zertifikat)), dateiname=a.zertifikat)
    except Exception as exc:
        stem = Path(a.zertifikat).stem
        r = {"discriminant": -int(stem) if stem.isdigit() else stem, "file": Path(a.zertifikat).name,
             "verified": False, "stage": "exception", "detail": repr(exc)}
    print(json.dumps(r, indent=1)); sys.exit(0 if r["verified"] else 1)


if __name__ == "__main__":
    main()
