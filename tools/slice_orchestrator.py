#!/usr/bin/env python3
"""Ein Einstieg fuer alle Slice-Rechnungen (p = 3, 5, 7).

Konsolidiert die Orchestrierung der drei historisch getrennten
Ketten: Stufenregister je Primzahl, artefaktbasiertes Resume,
Prozesspool ueber Koerper, Preflight-Checks und eine Bilanz, die
ausschliesslich aus dem Baum gezaehlt wird (nie aus Protokolldateien).
Die mathematischen Werkzeuge selbst (Builder, Verifier, Suchen)
bleiben unangetastet; dieses Skript ruft sie nur auf.

    python3 tools/slice_orchestrator.py --prime 5 \
        --fields experiments/p5-block23/fields.tsv \
        --tree   experiments/p5-block23/tree [--workers 4]
        [--stages cm,chars] [--only -541579031] [--dry-run]
        [--preflight-only]

fields.tsv-Schema (Tab-getrennt, Kopfzeile):
  p=5:  D_K  class_group  class_number  base_polynomial  radicand
  p=7:  D_K  class_group  base_polynomial
  p=3:  D_K  class_group  [weitere Spalten frei]

Design: experiments/pipeline-konsolidierung/DESIGN.md
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------------ Felder

@dataclass(frozen=True)
class Field:
    D: int
    class_group: str
    pol: str
    radicand: int
    extra: dict = field(default_factory=dict, compare=False)

    @property
    def absd(self) -> int:
        return abs(self.D)


def read_fields(path: Path, prime: int) -> list[Field]:
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        raise SystemExit(f"leere Koerperliste: {path}")
    out = []
    for r in rows:
        D = int(r["D_K"])
        pol = (r.get("base_polynomial") or "").strip()
        if not pol:
            a = (1 - D) // 4 if D % 4 == 1 else -D // 4
            pol = f"s^2-s+{a}" if D % 4 == 1 else f"s^2+{a}"
        rad = int(r.get("radicand") or (abs(D) if D % 4 == 1 else abs(D) // 4))
        out.append(Field(D, r.get("class_group", "").strip("[]"),
                         pol, rad, dict(r)))
    out.sort(key=lambda f: f.absd)
    return out


# ------------------------------------------------------------------ Stufen

@dataclass(frozen=True)
class Stage:
    name: str
    # done(tree, f) -> Path des Fertig-Artefakts (existiert+Muster = skip)
    artifact: object
    success_pattern: str | None
    run: object          # run(ctx, f) -> None (wirft bei Fehler)
    per_field: bool = True


@dataclass
class Ctx:
    prime: int
    tree: Path
    fields_path: Path
    env: dict
    log_dir: Path

    def field_dir(self, f: Field) -> Path:
        return self.tree / "fields" / f"D-{f.absd}"

    def log(self, stage: str, f: Field | None = None) -> Path:
        tag = f"{stage}-D{f.absd}" if f else stage
        return self.log_dir / f"{tag}.log"


def sh(ctx: Ctx, args, log: Path, cwd: Path | None = None,
       env_extra: dict | None = None, timeout: int | None = None):
    env = dict(os.environ)
    env.update(ctx.env)
    if env_extra:
        env.update(env_extra)
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as fh:
        fh.write(f"\n=== {time.strftime('%F %T')} $ {' '.join(map(str, args))}\n")
        fh.flush()
        rc = subprocess.run([str(a) for a in args], stdout=fh,
                            stderr=subprocess.STDOUT, cwd=cwd, env=env,
                            timeout=timeout).returncode
    if rc != 0:
        raise RuntimeError(f"exit {rc}: {' '.join(map(str, args))} (Log: {log})")


def artifact_ok(path: Path, pattern: str | None) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if pattern is None:
        return True
    try:
        return re.search(pattern, path.read_text(errors="replace")) is not None
    except OSError:
        return False


# ------------------------------------------------------- Preflight-Checks

@dataclass(frozen=True)
class Tool:
    label: str
    path: Path                 # Binary oder Skript
    sources: tuple = ()        # Quelldateien: Binary darf nicht aelter sein
    rebuild_hint: str = ""


def preflight(tools: list[Tool], fields: list[Field], tree: Path) -> list[str]:
    problems = []
    for t in tools:
        if not t.path.exists():
            problems.append(f"fehlt: {t.path}"
                            + (f"  ->  {t.rebuild_hint}" if t.rebuild_hint else ""))
            continue
        if t.sources:
            bin_m = t.path.stat().st_mtime
            newer = [s for s in t.sources
                     if Path(s).exists() and Path(s).stat().st_mtime > bin_m]
            if newer:
                problems.append(
                    f"veraltet: {t.path} ist aelter als "
                    + ", ".join(str(n) for n in newer)
                    + (f"  ->  {t.rebuild_hint}" if t.rebuild_hint else ""))
    for exe in ("gp",):
        if shutil.which(exe) is None and not (Path.home() / ".local/bin" / exe).exists():
            problems.append(f"nicht im PATH: {exe}")
    seen = set()
    for f in fields:
        if f.D in seen:
            problems.append(f"doppelte Diskriminante in fields.tsv: {f.D}")
        seen.add(f.D)
    try:
        tree.mkdir(parents=True, exist_ok=True)
        probe = tree / ".write-probe"
        probe.write_text("ok")
        probe.unlink()
    except OSError as e:
        problems.append(f"Baum nicht beschreibbar: {tree} ({e})")
    return problems


def gp_bin() -> str:
    local = Path.home() / ".local/bin/gp"
    return str(local) if local.exists() else "gp"


# ------------------------------------------------------------- p=5-Kette

def p5_stages(ctx: Ctx):
    cmdir = ROOT / "cm-constructor"
    verifier = ROOT / "certificate/verify_certificate"

    def out(f: Field) -> Path:
        return ctx.field_dir(f)

    def st_cm(ctx, f):
        d = out(f); d.mkdir(parents=True, exist_ok=True)
        sh(ctx, [cmdir / "cm_construct", f.pol, d / "cm.gp",
                 ctx.env.get("SAFETY_BITS", "768")], ctx.log("cm", f))

    def st_audit(ctx, f):
        d = out(f)
        sh(ctx, [cmdir / "audit_cm_fields", d / "cm.gp"], ctx.log("audit", f))
        txt = ctx.log("audit", f).read_text(errors="replace")
        if "6/6 VERIFIED" not in txt:
            raise RuntimeError(f"audit_cm_fields ohne 6/6 VERIFIED (D={f.D})")
        (d / "audit-ok").write_text("6/6 VERIFIED\n")

    def st_chars(ctx, f):
        d = out(f)
        for i in range(1, 7):
            cert = d / f"char{i}-cert.gp"
            mat = d / f"char{i}-mat.gp"
            if cert.is_file() and cert.stat().st_size and mat.is_file():
                continue
            sh(ctx, [cmdir / "cm_character_driver", "5", f.pol, str(i), mat],
               ctx.log(f"char{i}", f),
               env_extra={"MASSEY_CM_OUTPUT": str(d / "cm.gp"),
                          "MASSEY_CM_CHARACTER": str(i),
                          "MASSEY_CERTIFICATE_EXPORT": str(cert)})

    def st_assemble(ctx, f):
        d = out(f)
        sh(ctx, [sys.executable,
                 ROOT / "experiments/p5-slice2/assemble_certificate.py",
                 d / "certificate.gp",
                 *[d / f"char{i}-cert.gp" for i in range(1, 7)]],
           ctx.log("assemble", f))

    def st_finish(ctx, f):
        d = out(f)
        sh(ctx, [ROOT / "parallelization/finish_driver", "5", f.pol,
                 d / "result.gp", "exhaustive",
                 *[d / f"char{i}-mat.gp" for i in range(1, 7)]],
           ctx.log("finish", f))

    def st_verify(ctx, f):
        d = out(f)
        sh(ctx, [verifier, d / "certificate.gp", d / "result.gp"],
           ctx.log("verify", f))
        txt = ctx.log("verify", f).read_text(errors="replace")
        if "CERTIFICATE VERIFIED" not in txt:
            raise RuntimeError(f"Verifier ohne VERIFIED (D={f.D})")
        (d / "verify.txt").write_text(txt)

    tools = [
        Tool("cm_construct", cmdir / "cm_construct",
             (cmdir / "cm_construct.c",),
             f"make -C {cmdir} PARI=$HOME/.local"),
        Tool("cm_character_driver", cmdir / "cm_character_driver",
             tuple((ROOT / "src").glob("*.c"))
             + (ROOT / "parallelization/character_driver.c",),
             f"make PARI=$HOME/.local && make -C {cmdir} PARI=$HOME/.local"),
        Tool("finish_driver", ROOT / "parallelization/finish_driver",
             tuple((ROOT / "src").glob("*.c")),
             "make PARI=$HOME/.local && make -C parallelization PARI=$HOME/.local"),
        Tool("verify_certificate", verifier,
             (ROOT / "certificate/verify_certificate.c",),
             "make -C certificate PARI=$HOME/.local"),
    ]
    stages = [
        Stage("cm", lambda c, f: c.field_dir(f) / "cm.gp", None, st_cm),
        Stage("audit", lambda c, f: c.field_dir(f) / "audit-ok", None, st_audit),
        Stage("chars", lambda c, f: c.field_dir(f) / "char6-cert.gp", None, st_chars),
        Stage("assemble", lambda c, f: c.field_dir(f) / "certificate.gp", None, st_assemble),
        Stage("finish", lambda c, f: c.field_dir(f) / "result.gp", None, st_finish),
        Stage("verify", lambda c, f: c.field_dir(f) / "verify.txt",
              "CERTIFICATE VERIFIED", st_verify),
    ]
    return stages, tools


# ------------------------------------------------------------- p=7-Kette

def p7_stages(ctx: Ctx):
    sonde = ROOT / "experiments/p7-sondierung"
    verifier = ROOT / "examples/p3/verifier/verify_certificate"

    def st_cm(ctx, f):
        d = ctx.field_dir(f); d.mkdir(parents=True, exist_ok=True)
        sh(ctx, [sonde / "cm7/cm_construct7", f.pol, d / "cm7-fields-all.gp",
                 ctx.env.get("SAFETY_BITS", "768")], ctx.log("cm7", f))

    def st_build(ctx, f):
        d = ctx.field_dir(f)
        (d / "build-result").mkdir(exist_ok=True)
        sh(ctx, [gp_bin(), "-q", sonde / "builder7/build_certificate7.gp"],
           ctx.log("build7", f),
           env_extra={"P7_CM_FIELDS": str(d / "cm7-fields-all.gp"),
                      "P7_DISC": str(f.D),
                      "P7_EXPECTED_CYC": f"[{f.class_group}]",
                      "P7_RESULT_DIR": str(d / "build-result"),
                      "P7_CERT_PATH": str(d / "build-result/certificate.gp")})

    def st_transverse(ctx, f):
        d = ctx.field_dir(f)
        sh(ctx, [sys.executable, sonde / "transverse7.py",
                 d / "build-result/matrices.tsv"], ctx.log("transverse7", f))
        shutil.copyfile(ctx.log("transverse7", f), d / "transverse7.txt")

    def st_verify(ctx, f):
        d = ctx.field_dir(f)
        sh(ctx, [verifier, d / "build-result/certificate.gp"],
           ctx.log("verify", f))
        txt = ctx.log("verify", f).read_text(errors="replace")
        if "CERTIFICATE VERIFIED" not in txt:
            raise RuntimeError(f"Verifier ohne VERIFIED (D={f.D})")
        (d / "verify.txt").write_text(txt)

    tools = [
        Tool("cm_construct7", sonde / "cm7/cm_construct7",
             (sonde / "cm7/cm_construct_pdeg.c",),
             f"make -C {sonde}/cm7 PARI=$HOME/.local cm_construct7"),
        Tool("builder7", sonde / "builder7/build_certificate7.gp"),
        Tool("transverse7", sonde / "transverse7.py"),
        Tool("verify_certificate (generalisiert)", verifier,
             (ROOT / "examples/p3/verifier/verify_certificate.c",),
             "make -C examples/p3/verifier PARI=$HOME/.local"),
    ]
    stages = [
        Stage("cm7", lambda c, f: c.field_dir(f) / "cm7-fields-all.gp", None, st_cm),
        Stage("build7", lambda c, f: c.field_dir(f) / "build-result/certificate.gp",
              None, st_build),
        Stage("transverse7", lambda c, f: c.field_dir(f) / "transverse7.txt",
              "ERGEBNIS", st_transverse),
        Stage("verify", lambda c, f: c.field_dir(f) / "verify.txt",
              "CERTIFICATE VERIFIED", st_verify),
    ]
    return stages, tools


# ------------------------------------------------------------- p=3-Kette

def p3_stages(ctx: Ctx):
    """Phase 1: Adapter auf die vorhandenen p3-slice2-Runner (die selbst
    schon resumierbar und parametrisiert sind).  Der Orchestrator
    liefert Preflight und einheitliche Aufrufe; Feinintegration der
    Einzelstufen folgt in Phase 2."""
    slice2 = ROOT / "experiments/p3-slice2"

    def st_all(ctx, f=None):
        sh(ctx, ["sh", slice2 / "run_all.sh"], ctx.log("run_all"),
           cwd=slice2)

    tools = [
        Tool("p3-Pilot", ROOT / "experiments/p3-arithmetic-pilot/arithmetic.gp"),
        Tool("p3-Builder", ROOT / "examples/p3/builder/build_certificates.py"),
        Tool("p3-Verifier", ROOT / "examples/p3/verifier/verify_certificate",
             (ROOT / "examples/p3/verifier/verify_certificate.c",),
             "make -C examples/p3/verifier PARI=$HOME/.local"),
        Tool("block_sweep", ROOT / "examples/p3/strong-freeness/block_sweep.py"),
    ]
    stages = [Stage("run_all", lambda c, f: c.tree / "logs/run_all.log",
                    "FERTIG", st_all, per_field=False)]
    return stages, tools


# ------------------------------------------------------------ Ausfuehrung

def run_field(ctx: Ctx, stages, f: Field, wanted: set[str] | None):
    done, ran = [], []
    for st in stages:
        if not st.per_field:
            continue
        if wanted and st.name not in wanted:
            continue
        art = st.artifact(ctx, f)
        if artifact_ok(art, st.success_pattern):
            done.append(st.name)
            continue
        st.run(ctx, f)
        if not artifact_ok(art, st.success_pattern):
            raise RuntimeError(
                f"Stufe {st.name}: Artefakt fehlt/unvollstaendig: {art}")
        ran.append(st.name)
    return f.D, done, ran


def summarize(ctx: Ctx, stages, fields) -> str:
    lines = [f"Bilanz (aus dem Baum gezaehlt, {time.strftime('%F %T')}):"]
    for st in stages:
        if not st.per_field:
            continue
        n = sum(1 for f in fields
                if artifact_ok(st.artifact(ctx, f), st.success_pattern))
        lines.append(f"  {st.name:12s} {n}/{len(fields)}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", type=int, required=True, choices=(3, 5, 7))
    ap.add_argument("--fields", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--stages", type=str, default=None,
                    help="Kommagetrennte Stufen; Default: alle")
    ap.add_argument("--only", type=int, nargs="*", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--preflight-only", action="store_true")
    a = ap.parse_args(argv)

    fields = read_fields(a.fields, a.prime)
    if a.only:
        wanted_d = {abs(d) for d in a.only}
        fields = [f for f in fields if f.absd in wanted_d]
        if not fields:
            raise SystemExit("--only trifft keinen Koerper der Liste")

    tree = a.tree.resolve()
    ctx = Ctx(a.prime, tree, a.fields.resolve(),
              {k: v for k, v in os.environ.items()
               if k.startswith(("MASSEY_", "P7_", "P3_", "SAFETY_"))},
              tree / "logs")
    ctx.log_dir.mkdir(parents=True, exist_ok=True)

    stages, tools = {3: p3_stages, 5: p5_stages, 7: p7_stages}[a.prime](ctx)
    wanted = set(a.stages.split(",")) if a.stages else None
    if wanted:
        unknown = wanted - {s.name for s in stages}
        if unknown:
            raise SystemExit(f"unbekannte Stufen: {sorted(unknown)}")

    problems = preflight(tools, fields, tree)
    if problems:
        print("PREFLIGHT-PROBLEME:")
        for p in problems:
            print("  -", p)
        return 2
    print(f"Preflight ok: {len(fields)} Koerper, Stufen "
          f"{[s.name for s in stages]}, Baum {tree}")
    if a.preflight_only:
        return 0
    if a.dry_run:
        for f in fields:
            todo = [s.name for s in stages if s.per_field
                    and (not wanted or s.name in wanted)
                    and not artifact_ok(s.artifact(ctx, f), s.success_pattern)]
            print(f"  D={f.D}: " + (", ".join(todo) if todo else "fertig"))
        return 0

    global_stage = [s for s in stages if not s.per_field]
    if global_stage:
        for st in global_stage:
            art = st.artifact(ctx, None)
            if artifact_ok(art, st.success_pattern):
                print(f"Stufe {st.name}: Artefakt vorhanden, uebersprungen")
            else:
                st.run(ctx)
        print(summarize(ctx, stages, fields))
        return 0

    failures = []
    started = time.monotonic()
    if a.workers > 1:
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(run_field, ctx, stages, f, wanted): f
                    for f in fields}
            for fut, f in futs.items():
                try:
                    D, done, ran = fut.result()
                    if ran:
                        print(f"D={D}: neu {','.join(ran)}", flush=True)
                except Exception as e:
                    failures.append((f.D, str(e)))
                    print(f"D={f.D}: FEHLER {e}", flush=True)
    else:
        for f in fields:
            try:
                D, done, ran = run_field(ctx, stages, f, wanted)
                if ran:
                    print(f"D={D}: neu {','.join(ran)}", flush=True)
            except Exception as e:
                failures.append((f.D, str(e)))
                print(f"D={f.D}: FEHLER {e}", flush=True)

    print(summarize(ctx, stages, fields))
    print(f"Dauer {time.monotonic()-started:.0f}s; Fehler: {len(failures)}")
    for D, msg in failures:
        print(f"  D={D}: {msg}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
