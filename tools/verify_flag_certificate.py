#!/usr/bin/env python3
"""Standalone verifier for flag certificates (plain Python, no computer algebra system).

For a certificate (records/p3/flag-certificates/ or records/p5/flag-certificates/) it checks
  1. that the words avoiding the head-word set W are counted by
     1/(1-3z+3z^3) in every degree (prefix automaton, Lemma A.4);
  2. that after the change of variables over the stated field F_q,
     q = p^degree, the completion of the three cubic relations
     (Lemma A.5) has exactly the head words W and every overlap up to
     degree 2*max|W|-1 reduces to zero (finite completed basis).
Together this is strong freeness (Corollary A.2), hence mildness
(Lemma 3.7, Proposition 3.8).  The relation tensor of the field is read
from the verified records of the repository (p = 3: tensor_3_by_27 of
the verification records; p = 5: cubic_relation_matrix of result.gp).

  python3 verify_flag_certificate.py <certificate.json> [--details]
  python3 verify_flag_certificate.py --all <directory> [--jobs N] [--out report.jsonl] [--no-resume]
"""
import argparse, json, os, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flag_cell_checker as zelle  # noqa: E402

# Tensor sources: the verified census records of the repository (key tensor_3_by_27)
REPO = HERE.parent
STANDARD_TENSOREN = [REPO / "records" / "p3" / "results" / name for name in (
    "verification-004.json", "verification-005.json",
    "slice2/verification-starthinker.json", "block23-verification-starthinker.json")]
P5_RESULT = REPO / "records" / "p5" / "D-{absD}" / "result.gp"   # cubic_relation_matrix over F_5
_TENSOR_CACHE = {}

def _lade_result_gp(p, prime, D):
    import re
    text = p.read_text()
    mp = re.search(r'"p",\s*(\d+)', text)
    md = re.search(r'"base_discriminant",\s*(-\d+)', text)
    if not mp or int(mp.group(1)) != prime: raise KeyError(f"{p}: prime of the record differs from the certificate")
    if not md or int(md.group(1)) != D: raise KeyError(f"{p}: base_discriminant of the record differs from the certificate")
    m = re.search(r'"cubic_relation_matrix",\s*\[(.*?)\]\]', text, re.S)
    if not m: raise KeyError(f"{p}: no cubic_relation_matrix")
    rows = [r for r in m.group(1).replace("[", "").replace("]", "").split(";") if r.strip()]
    T = [[int(v) % prime for v in r.split(",")] for r in rows]
    if len(T) != 3 or any(len(r) != 27 for r in T): raise KeyError(f"{p}: tensor is not 3 x 27")
    return T

def _lade_records(p):
    data = json.loads(p.read_text())
    tens = {}
    for r in data.get("records", []):
        if r.get("status") == "VERIFIED" and "tensor_3_by_27" in r:
            tens[str(r["discriminant"])] = r["tensor_3_by_27"]
    return tens

def lade_tensor(D, pfade, prime=3):
    key = str(D)
    if prime == 5:
        p = Path(str(P5_RESULT).format(absD=abs(D)))
        if not p.exists(): raise KeyError(f"tensor of {D}: {p} not found")
        return _lade_result_gp(p, 5, D), str(p)
    for p in pfade:
        p = Path(p)
        if not p.exists(): continue
        if p not in _TENSOR_CACHE:
            _TENSOR_CACHE[p] = _lade_records(p) if "records" in p.parts else json.loads(p.read_text())
        if key in _TENSOR_CACHE[p]: return _TENSOR_CACHE[p][key], str(p)
    raise KeyError(f"tensor of {D} not found in the census records")

VERIFIER_VERSION = "2026-09-01.2"
MODULI = {(prime, e): list(m) for (prime, e), m in zelle.FiniteField.MODULI.items()}
PRIMZAHLEN = (3, 5)
PFLICHT = {"discriminant", "field", "basis", "head_words", "proper_prefixes", "search", "claim", "alternatives"}

def schema(cert, dateiname=None):
    """Strict format check; returns (prime, degree, canonical basis) or raises ValueError."""
    if set(cert) != PFLICHT:
        raise ValueError(f"keys: expected {sorted(PFLICHT)}, found {sorted(cert)}")
    D = cert["discriminant"]
    if not isinstance(D, int) or D >= 0: raise ValueError("discriminant must be a negative integer")
    if dateiname is not None and Path(dateiname).name != f"{abs(D)}.json":
        raise ValueError("file name does not match the discriminant")
    fld = cert["field"]
    if set(fld) != {"q", "degree", "modulus_low_to_high"}: raise ValueError("field keys")
    e = fld["degree"]; q = fld["q"]
    if not (isinstance(e, int) and e >= 1 and isinstance(q, int)): raise ValueError("degree and q must be positive integers")
    prime = next((pr for pr in PRIMZAHLEN if pr ** e == q), None)
    if prime is None: raise ValueError("q must equal p^degree for p in {3, 5}")
    if (prime, e) not in MODULI: raise ValueError(f"no documented modulus for p={prime}, degree {e}")
    if list(fld["modulus_low_to_high"]) != MODULI[(prime, e)]: raise ValueError("modulus differs from the documented one")
    B = cert["basis"]
    if not (isinstance(B, list) and len(B) == 3 and all(isinstance(r, list) and len(r) == 3 for r in B)):
        raise ValueError("basis must be a 3x3 matrix")
    kan = []
    for row in B:
        krow = []
        for x in row:
            if e == 1:
                if not (isinstance(x, int) and 0 <= x < prime): raise ValueError(f"basis entries over F_{prime} must be 0..{prime-1}")
                krow.append([x])
            else:
                if not (isinstance(x, list) and len(x) == e and all(isinstance(c, int) and 0 <= c < prime for c in x)):
                    raise ValueError(f"basis entries must be coefficient lists of length {e} with entries 0..{prime-1}")
                krow.append(list(x))
        kan.append(krow)
    W = cert["head_words"]
    if not (isinstance(W, list) and W and all(isinstance(w, str) and w and set(w) <= set("123") for w in W)):
        raise ValueError("head_words must be nonempty words over 1,2,3")
    if len(set(W)) != len(W): raise ValueError("duplicate head words")
    if any(u != w and u in w for u in W for w in W): raise ValueError("head words must be factor-free")
    s_erw = len({""} | {w[:i] for w in W for i in range(1, len(w))})
    if cert["proper_prefixes"] != s_erw: raise ValueError(f"proper_prefixes must be {s_erw}")
    if not isinstance(cert["claim"], str): raise ValueError("claim must be a string (it is documentation; the verifier does not interpret it)")
    return prime, e, kan

def determinante_nichtnull(kan, field):
    M = [[field.elt(x) for x in row] for row in kan]
    def m2(a, b, c, d): return field.sub(field.mul(a, d), field.mul(b, c))
    det = field.add(field.add(field.mul(M[0][0], m2(M[1][1], M[1][2], M[2][1], M[2][2])),
                              field.neg(field.mul(M[0][1], m2(M[1][0], M[1][2], M[2][0], M[2][2])))),
                    field.mul(M[0][2], m2(M[1][0], M[1][1], M[2][0], M[2][1])))
    return not field.is_zero(det)

def pruefe(cert, pfade, details=False, dateiname=None):
    D = cert["discriminant"]
    try:
        prime, e, kan = schema(cert, dateiname)
    except ValueError as exc:
        return {"discriminant": D, "verified": False, "stage": "schema", "detail": str(exc)}
    field = zelle.FiniteField(e, prime)
    if not determinante_nichtnull(kan, field):
        return {"discriminant": D, "prime": prime, "field": prime**e, "verified": False, "stage": "basis-singular"}
    tensor, quelle = lade_tensor(D, pfade, prime)
    t0 = time.monotonic()
    transformed = zelle.transform_tensor(tensor, kan, field)
    res = zelle.test_transformed(transformed, field, cert["head_words"], details=details, check_termination=True)
    ok = bool(res.get("W_regular", False))
    out = {"discriminant": D, "prime": prime, "field": prime**e, "verified": ok, "verifier_version": VERIFIER_VERSION,
           "head_words": len(cert["head_words"]), "max_head_degree": max(map(len, cert["head_words"])),
           "milliseconds": round(1000*(time.monotonic()-t0), 1), "tensor_source": quelle}
    if not ok: out["stage"] = res.get("failed_stage") or res.get("failed"); out["detail"] = {k: v for k, v in res.items() if k not in ("trace",)}
    return out

def digest(pfad, pfade):
    """Digest aus Zertifikatsdatei, Tensor und Verifier-Version (fuer das Resume)."""
    import hashlib
    h = hashlib.sha256(Path(pfad).read_bytes())
    try:
        cert = json.load(open(pfad)); fld = cert["field"]
        prime = next((pr for pr in PRIMZAHLEN if pr ** fld["degree"] == fld["q"]), 3)
        t, _ = lade_tensor(cert["discriminant"], pfade, prime)
        h.update(json.dumps(t, separators=(",", ":")).encode())
    except Exception:
        h.update(b"no-tensor")
    h.update(VERIFIER_VERSION.encode())
    return h.hexdigest()

def _job(arg):
    pfad, pfade = arg
    try:
        r = pruefe(json.load(open(pfad)), pfade, dateiname=pfad)
    except Exception as exc:
        r = {"discriminant": Path(pfad).stem, "verified": False, "stage": "exception", "detail": repr(exc)}
    r["file"] = Path(pfad).name; r["digest"] = digest(pfad, pfade)
    return r

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("zertifikat", nargs="?")
    ap.add_argument("--all"); ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", default="flag-certificates-verification.jsonl")
    ap.add_argument("--details", action="store_true")
    ap.add_argument("--no-resume", action="store_true", help="recheck every certificate, ignore the report")
    a = ap.parse_args()
    pfade = list(STANDARD_TENSOREN)
    if a.zertifikat:
        try:
            r = pruefe(json.load(open(a.zertifikat)), pfade, a.details, dateiname=a.zertifikat)
        except Exception as exc:
            stem = Path(a.zertifikat).stem
            r = {"discriminant": -int(stem) if stem.isdigit() else stem, "file": Path(a.zertifikat).name,
                 "verified": False, "stage": "exception", "detail": repr(exc)}
        print(json.dumps(r, indent=1)); sys.exit(0 if r["verified"] else 1)
    dateien = sorted(str(p) for p in Path(a.all).glob("*.json"))
    # Resume nur fuer echte Erfolge mit identischem Digest (Zertifikat + Tensor + Verifier-Version)
    erledigt = {}
    if os.path.exists(a.out) and not a.no_resume:
        for line in open(a.out):
            try:
                r = json.loads(line)
                if r.get("verified") and r.get("digest"): erledigt[r["file"]] = r
            except Exception: pass
    offen, uebernommen = [], []
    for f in dateien:
        alt = erledigt.get(Path(f).name)
        if alt and alt["digest"] == digest(f, pfade): uebernommen.append(alt)
        else: offen.append(f)
    print(f"{len(dateien)} certificates, {len(uebernommen)} unchanged and already verified, {len(offen)} to check, {a.jobs} jobs", flush=True)
    ergebnisse = list(uebernommen)
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, r in enumerate(ex.map(_job, [(f, pfade) for f in offen], chunksize=8), 1):
            ergebnisse.append(r)
            if not r["verified"]: print("  !! NOT VERIFIED:", r["file"], r.get("stage"), r.get("detail", ""), flush=True)
            if i % 500 == 0: print(f"  {i}/{len(offen)}", flush=True)
    # Gesamtbericht ueber ALLE aktuellen Dateien neu schreiben
    ergebnisse.sort(key=lambda r: abs(int(r["discriminant"])) if str(r["discriminant"]).lstrip("-").isdigit() else 0)
    with open(a.out, "w") as fh:
        for r in ergebnisse: fh.write(json.dumps(r) + "\n")
    ok = sum(1 for r in ergebnisse if r["verified"]); print(f"total: {len(ergebnisse)} certificates, {ok} verified, {len(ergebnisse)-ok} not verified", flush=True)
    sys.exit(0 if ok == len(ergebnisse) == len(dateien) else 1)

if __name__ == "__main__":
    main()
