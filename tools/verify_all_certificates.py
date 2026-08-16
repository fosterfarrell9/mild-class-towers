#!/usr/bin/env python3
"""Verify committed certificates and report the outcome per field.

Certificates are meant to be checkable anywhere, not only where they were
written -- that is what recording their integral bases is for.  Running this
on a machine that produced none of them is what confirms it, and keeps
confirming it.

Run from the repository root, with the verifier built:

    make -C verifier PARI=$HOME/.local
    python3 tools/verify_all_certificates.py

With no options the p = 5 collection is checked field by field, each
certificate cross-checked against its committed result record, one line of
output per field.  Two options change that:

    --all        check the p = 3, p = 5 and p = 7 collections -- every
                 certificate.gp found under each -- and report counts per
                 prime instead of a line per field.
    --read-only  give the verifier nothing but the certificate itself, and
                 its committed factorization hints where they exist; run it
                 from a scratch directory outside the repository.

Either option selects the counting sweep, and --all only widens that sweep
from p = 5 to all three collections.  The full check is

    python3 tools/verify_all_certificates.py --all --read-only

Nothing here writes to the repository under any option: this harness opens
no file for writing, and the verifier opens its inputs for reading only.
What --read-only adds is that the certificate is then judged alone, with no
result record beside it, from a working directory where a stray relative
path could not reach the tree in the first place.

Each certificate is given a bounded number of seconds; a run that reaches
the limit is reported as such rather than left to hang, since a misread
basis can send the ideal arithmetic into an integer factorization with no
useful bound.  The exit status is nonzero unless every certificate verified.
"""

import argparse
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TIMEOUT = "300"

PRIMES = (3, 5, 7)

# The p = 7 certificates run orders of magnitude longer than the p = 3 ones,
# so the sweep cannot share the per-field limit above: measured on one 32-core
# machine, a p = 3 certificate takes about two seconds and a p = 5 one about
# twenty-five, while K-931506071-p7 takes half an hour.  The limit is set well
# above that worst case rather than just above it, so that a slower machine
# reports what it found instead of a timeout it only produced by being slow.
SWEEP_TIMEOUT = 7200

PROGRESS_EVERY = 500


def result_records(root):
    """Map discriminant to result record, by content: the example directories
    are not named consistently enough to go by path."""
    records = {}
    for path in root.glob("records/p5/**/result.gp"):
        match = re.search(
            r'"base_discriminant", *(-?\d+)', path.read_text()[:3000])
        if match:
            records[match.group(1)] = path
    return records


def certificates(collection):
    """Every certificate in a collection, at whatever depth it sits.

    The trees are not shaped alike: p = 3 shards its fields into block
    directories, p = 5 and p = 7 keep one directory per field directly under
    the collection.  Searching for the file rather than for a fixed depth
    means a re-sharded collection is still swept in full."""
    return sorted(collection.rglob("certificate.gp"))


def check(certificate, verifier, workdir, timeout):
    """Run the verifier on one certificate and classify what came back.

    No result record is passed.  The point of this sweep is the claim that a
    certificate stands on its own, and passing a record would let a run pass
    on agreement with something the certificate did not say.  Committed
    factorization hints are passed where they exist: the verifier proves
    every hint prime before adding it, so a hint can only make the run
    faster or make it stop, never make it accept more."""
    command = [str(verifier), str(certificate)]
    hints = certificate.parent / "hints.gp"
    if hints.exists():
        command.append(str(hints))
    try:
        proc = subprocess.run(command, cwd=workdir, capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    if proc.returncode != 0:
        return "FAILED"
    output = proc.stdout + proc.stderr
    return "VERIFIED" if "CERTIFICATE VERIFIED" in output else "FAILED"


def sweep(root, verifier, primes, timeout, limit):
    """Check whole collections and report how many of each verified."""
    counts = []
    failures = []
    # A scratch directory in the system temporary area, never the tree: it
    # makes writing nothing into the repository structural rather than a
    # promise, since a relative path opened by anything below cannot reach
    # the tree from here.  Nothing is put in it; it exists to be somewhere
    # else.
    with tempfile.TemporaryDirectory(prefix="verify-certificates-") as workdir:
        for prime in primes:
            collection = root / "certificates" / f"p{prime}"
            if not collection.is_dir():
                sys.exit(f"no such collection: {collection}")
            found = certificates(collection)
            if limit is not None:
                found = found[:limit]
            print(f"p{prime}: {len(found)} certificate(s) under "
                  f"{collection.relative_to(root)}", flush=True)

            started = time.monotonic()
            passed = 0
            for number, certificate in enumerate(found, start=1):
                status = check(certificate, verifier, workdir, timeout)
                if status == "VERIFIED":
                    passed += 1
                else:
                    path = certificate.relative_to(root)
                    failures.append((status, path))
                    print(f"  {status}: {path}", flush=True)
                if number % PROGRESS_EVERY == 0 or number == len(found):
                    elapsed = time.monotonic() - started
                    print(f"  {number}/{len(found)} checked, "
                          f"{number - passed} failed, {elapsed:.0f}s",
                          flush=True)
            counts.append((prime, len(found), passed))

    print(f"\n{'collection':<12}{'checked':>9}{'passed':>9}{'failed':>9}")
    for prime, checked, passed in counts:
        print(f"{f'p{prime}':<12}{checked:>9}{passed:>9}{checked - passed:>9}")
    checked = sum(row[1] for row in counts)
    passed = sum(row[2] for row in counts)
    print(f"{'total':<12}{checked:>9}{passed:>9}{checked - passed:>9}")

    if failures:
        print(f"\n{len(failures)} certificate(s) did not verify:",
              file=sys.stderr)
        for status, path in failures:
            print(f"  {status:<9}{path}", file=sys.stderr)
    return 1 if failures else 0


def verify_p5_fields(root, verifier):
    """Check the p = 5 collection field by field against its result records."""
    cert_dir = root / "certificates" / "p5"
    records = result_records(root)

    failures = 0
    print(f"{'certificate':<16}{'status':<18}{'seconds':>8}"
          f"{'entries':>9}  cross-checked against")
    for directory in sorted(cert_dir.glob("K-*-p5")):
        certificate = directory / "certificate.gp"
        if not certificate.exists():
            continue
        head = certificate.open().read(400)
        match = re.match(r"\[\d+,\d+,\d+,[^,]+,(-?\d+),", head)
        if not match:
            print(f"{directory.name:<16}{'UNPARSED':<18}{0:>8}"
                  f"{'-':>9}  -", flush=True)
            failures += 1
            continue
        record = records.get(match.group(1))

        command = ["timeout", "-k", "10", TIMEOUT, str(verifier),
                   f"{directory.name}/certificate.gp"]
        if record:
            command.append(str(Path("..") / record.relative_to(root)))

        started = time.monotonic()
        proc = subprocess.run(command, cwd=cert_dir, capture_output=True,
                              text=True)
        elapsed = time.monotonic() - started
        output = proc.stdout + proc.stderr

        # How many entries there are is a property of the certificate -- the
        # principal example carries the nine doubled-character entries as
        # well -- so it is counted in the file, not in the output.  Counting
        # it in the output would make the denominator agree with the numerator
        # by construction and hide exactly the case worth seeing: a run that
        # stopped early.
        entries = output.count("AC1=PASS")
        expected = sum(1 for line in certificate.open()
                       if re.match(r"\[+\"", line))

        if proc.returncode in (124, 137):
            status = "TIMEOUT"
        elif "CERTIFICATE VERIFIED" not in output:
            status = "FAILED"
        elif record and "RESULT_RECORD_MATCH=PASS" not in output:
            status = "VERIFIED_NO_MATCH"
        else:
            status = "VERIFIED"
        if status != "VERIFIED":
            failures += 1

        against = record.parent.name if record else "nothing"
        print(f"{directory.name:<16}{status:<18}{elapsed:>8.1f}"
              f"{f'{entries}/{expected}':>9}  {against}", flush=True)

    if failures:
        print(f"\n{failures} certificate(s) did not verify", file=sys.stderr)
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(
        description="Verify committed certificates and report the outcome.")
    parser.add_argument(
        "--all", action="store_true",
        help="check the p = 3, p = 5 and p = 7 collections, not p = 5 alone")
    parser.add_argument(
        "--read-only", action="store_true",
        help="give the verifier the certificate alone, from a scratch "
             "directory outside the repository")
    parser.add_argument(
        "--timeout", type=int, default=None, metavar="SECONDS",
        help=f"per-certificate limit for the sweep (default {SWEEP_TIMEOUT})")
    parser.add_argument(
        "--limit", type=int, default=None, metavar="N",
        help="stop after the first N certificates of each collection; a "
             "smoke test, not a check of the collection")
    args = parser.parse_args()
    sweeping = args.all or args.read_only
    for option in ("timeout", "limit"):
        if getattr(args, option) is not None and not sweeping:
            parser.error(f"--{option} applies to --all/--read-only only")

    root = Path.cwd()
    verifier = root / "verifier" / "verify_certificate"
    if not verifier.exists():
        sys.exit("build the verifier first: "
                 "make -C verifier PARI=<pari-prefix>")
    if sweeping:
        timeout = SWEEP_TIMEOUT if args.timeout is None else args.timeout
        return sweep(root, verifier, PRIMES if args.all else (5,),
                     timeout, args.limit)
    return verify_p5_fields(root, verifier)


if __name__ == "__main__":
    sys.exit(main())
