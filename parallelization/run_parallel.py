#!/usr/bin/env python3
"""Per-character parallel computation of one field.

Spawns six character_driver processes (one prescribed character each, all
with the exact audit and inline certificate export), merges the six partial
certificates into one, assembles the result record with finish_driver, and
verifies the merged certificate with the shared standalone verifier against
the fresh result record (RESULT_RECORD_MATCH).

The proven sequential code in src/ is reused unmodified (the drivers link
against the main build's object files).

Usage, from this directory (repository built, verifier built):

  python3 run_parallel.py --polynomial 's^2-s+15260177' \
      --workdir work/D-61040707 [--limit exhaustive] \
      [--verify-against ../records/p5/batch-block0-01/D-61040707/result.gp]

--verify-against runs the regression comparison: the fresh result must
reproduce the committed record's verified entries, and the merged
certificate must also verify against the committed record.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "tools"))

from run_mildness_batch import result_entry, CERT_COMPARE_KEYS  # noqa: E402

LABELS = ["a", "b", "c", "a+b", "a+c", "b+c"]


def run_field(polynomial, workdir, limit, indices=None, finish=True,
              driver="character_driver"):
    workdir.mkdir(parents=True, exist_ok=True)
    procs = []
    started = time.monotonic()
    for i in (indices or range(1, 7)):
        env = dict(os.environ)
        env["MASSEY_CERTIFICATE_EXPORT"] = str(
            workdir / f"cert-{i}.gp")
        log = (workdir / f"char-{i}.log").open("w")
        child = subprocess.Popen(
            ["/usr/bin/time", "-v",
             str(HERE / driver), "5", polynomial, str(i),
             str(workdir / f"mat-{i}.gp")],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        stamper = subprocess.Popen(
            ["awk", '{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0;'
             ' fflush() }'],
            stdin=child.stdout, stdout=log)
        child.stdout.close()
        procs.append((i, child, stamper, log))
        print(f"spawned character {i} ({LABELS[i-1]})", flush=True)
    failures = []
    for i, proc, stamper, log in procs:
        rc = proc.wait()
        stamper.wait()
        log.close()
        peak = ""
        peak_kb = ""
        for line in (workdir / f"char-{i}.log").read_text().splitlines():
            if "Maximum resident set size" in line:
                peak_kb = int(line.rsplit(":", 1)[1])
                peak = f", peak RSS {peak_kb / 1048576:.1f} GiB"
        print(f"character {i} ({LABELS[i-1]}) exited {rc} "
              f"after {time.monotonic()-started:.0f}s{peak}", flush=True)
        record_stats(workdir, polynomial, i, rc,
                     round(time.monotonic() - started), peak_kb)
        if rc != 0:
            failures.append(i)
    if failures:
        raise SystemExit(f"character processes failed: {failures}")
    if not finish:
        print("subset complete; skipping merge/finish", flush=True)
        return

    merge_certificates(
        [workdir / f"cert-{i}.gp" for i in range(1, 7)],
        workdir / "certificate.gp")

    finish_log = (workdir / "finish.log").open("w")
    rc = subprocess.run(
        [str(HERE / "finish_driver"), "5", polynomial,
         str(workdir / "result.gp"), limit]
        + [str(workdir / f"mat-{i}.gp") for i in range(1, 7)],
        stdout=finish_log, stderr=subprocess.STDOUT).returncode
    finish_log.close()
    if rc != 0:
        raise SystemExit(f"finish_driver exited {rc}")
    print("result assembled", flush=True)


def record_stats(workdir, polynomial, index, rc, seconds, peak_kb):
    """Append one per-character cost record to the empirical table."""
    stats = HERE / "stats.tsv"
    if not stats.exists():
        stats.write_text("field\tpolynomial\tcharacter\tlabel\texit"
                         "\tseconds\tpeak_rss_kb\tstack_cap\n")
    cap = os.environ.get("MASSEY_PARISTACK_MAX", "default")
    with stats.open("a") as fh:
        fh.write(f"{workdir.name}\t{polynomial}\t{index}"
                 f"\t{LABELS[index-1]}\t{rc}\t{seconds}"
                 f"\t{peak_kb}\t{cap}\n")


def merge_certificates(paths, out):
    header = None
    entries = []
    for path in paths:
        text = path.read_text()
        idx = text.find(",\n[")
        assert idx > 0, f"{path}: malformed partial certificate"
        head = text[:idx]
        if header is None:
            header = head
        else:
            assert head == header, f"{path}: header differs"
        blob = text[idx + 3:]
        assert blob.endswith("\n]]\n"), f"{path}: missing finalizer"
        blob = blob[: -len("\n]]\n")]
        parts = blob.split(",\n")
        assert len(parts) == 3, f"{path}: expected three entries"
        entries.extend(parts)
    out.write_text(header + ",\n[" + ",\n".join(entries) + "\n]]\n")
    print(f"merged certificate: {out} ({len(entries)} entries)",
          flush=True)


def verify(certificate, result):
    proc = subprocess.run(
        [str(ROOT / "verifier" / "verify_certificate"),
         str(Path(certificate).resolve()), str(Path(result).resolve())],
        cwd=ROOT / "certificates" / "p5", capture_output=True, text=True)
    tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
    ok = (proc.returncode == 0
          and any("CERTIFICATE VERIFIED" in l for l in tail)
          and any("RESULT_RECORD_MATCH=PASS" in l for l in tail))
    return ok, tail


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--polynomial", required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--limit", default="exhaustive")
    parser.add_argument("--verify-against", type=Path, default=None)
    parser.add_argument("--characters", default=None,
                        help="comma list of indices 1..6: compute only "
                        "these and stop (no merge/finish)")
    parser.add_argument("--driver", default="character_driver",
                        help="driver binary: character_driver (default)"
                        " or character_driver_mt for the threaded build")
    parser.add_argument("--merge-only", action="store_true",
                        help="skip computation; merge existing partial "
                        "certificates, assemble and verify the result")
    args = parser.parse_args()

    if args.characters:
        indices = [int(v) for v in args.characters.split(",")]
        run_field(args.polynomial, args.workdir, args.limit,
                  indices=indices, finish=False, driver=args.driver)
        return
    if args.merge_only:
        merge_certificates(
            [args.workdir / f"cert-{i}.gp" for i in range(1, 7)],
            args.workdir / "certificate.gp")
        finish_log = (args.workdir / "finish.log").open("w")
        rc = subprocess.run(
            [str(HERE / "finish_driver"), "5", args.polynomial,
             str(args.workdir / "result.gp"), args.limit]
            + [str(args.workdir / f"mat-{i}.gp") for i in range(1, 7)],
            stdout=finish_log, stderr=subprocess.STDOUT).returncode
        finish_log.close()
        if rc != 0:
            raise SystemExit(f"finish_driver exited {rc}")
    else:
        run_field(args.polynomial, args.workdir, args.limit,
                  driver=args.driver)

    certificate = args.workdir / "certificate.gp"
    fresh = args.workdir / "result.gp"
    ok, tail = verify(certificate, fresh)
    print("verification vs fresh result:", tail)
    if not ok:
        raise SystemExit("verification against fresh result FAILED")

    if args.verify_against:
        committed = args.verify_against.resolve()
        old = committed.read_text()
        new = fresh.read_text()
        for key in CERT_COMPARE_KEYS:
            assert result_entry(old, key) == result_entry(new, key), \
                f"REGRESSION: {key} differs from committed record"
        print("fresh result matches committed record entry for entry")
        ok, tail = verify(certificate, committed)
        print("verification vs committed result:", tail)
        if not ok:
            raise SystemExit("verification against committed FAILED")
        print("REGRESSION PASS")

    print("PARALLEL FIELD COMPLETE")


if __name__ == "__main__":
    main()
