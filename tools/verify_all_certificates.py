#!/usr/bin/env python3
"""Verify every committed certificate and report one line per field.

Certificates are meant to be checkable anywhere, not only where they were
written -- that is what recording their integral bases is for.  Running this
on a machine that produced none of them is what confirms it, and keeps
confirming it.

Run from the repository root, with the verifier built:

    make -C verifier PARI=$HOME/.local
    python3 tools/verify_all_certificates.py

Each certificate is given at most TIMEOUT seconds; a run that reaches the
limit is reported as such rather than left to hang, since a misread basis can
send the ideal arithmetic into an integer factorization with no useful bound.
The exit status is nonzero unless every certificate verified.
"""

import re
import subprocess
import sys
import time
from pathlib import Path

TIMEOUT = "300"


def result_records(root):
    """Map discriminant to result record, by content: the example directories
    are not named consistently enough to go by path."""
    records = {}
    for path in root.glob("examples/p5/**/result.gp"):
        match = re.search(
            r'"base_discriminant", *(-?\d+)', path.read_text()[:3000])
        if match:
            records[match.group(1)] = path
    return records


def main():
    root = Path.cwd()
    cert_dir = root / "certificate"
    verifier = root / "verifier" / "verify_certificate"
    if not verifier.exists():
        sys.exit("build the verifier first: "
                 "make -C verifier PARI=<pari-prefix>")
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


if __name__ == "__main__":
    sys.exit(main())
