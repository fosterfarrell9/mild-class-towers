#!/usr/bin/env python3
"""Verify every p=3 certificate and persist reconstructed tensors.

The C verifier reads only the certificate.  This harness supplies timeouts,
never overwrites logs, and turns the verifier's tensor line into JSON for the
finite strong-freeness checks.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
VERIFIER = ROOT / "verifier" / "verify_certificate"
CERTIFICATES = ROOT / "certificates" / "p3"
RESULTS = ROOT / "records" / "p3" / "results"
TIMEOUT_SECONDS = 300

def bucket(discriminant: int) -> str:
    """Shard directory: floor(|D|/10^7), three digits."""
    return f"{abs(discriminant) // 10**7:03d}"


def source_tensor_path(discriminant: int) -> Path:
    """The independent baseline tensor for a field, by discriminant."""
    return (ROOT / "records" / "p3" / "source-tensors" / bucket(discriminant)
            / f"D-{abs(discriminant)}/tensor.json")


def next_log(directory: Path) -> tuple[int, Path]:
    number = 1
    while (directory / f"verify-{number:03d}.log").exists():
        number += 1
    return number, directory / f"verify-{number:03d}.log"


def next_result() -> Path:
    results = RESULTS
    results.mkdir(exist_ok=True)
    number = 1
    while (results / f"verification-{number:03d}.json").exists():
        number += 1
    return results / f"verification-{number:03d}.json"


def verify(directory: Path) -> dict[str, object]:
    certificate = directory / "certificate.gp"
    discriminant = -int(re.fullmatch(r"K-(\d+)-p3", directory.name).group(1))
    run, log_path = next_log(directory)
    started = time.monotonic()
    command = [str(VERIFIER), str(certificate)]
    hints = directory / "hints.gp"
    if hints.exists():
        command.append(str(hints))
    try:
        process = subprocess.run(
            command, cwd=HERE, text=True,
            capture_output=True, timeout=TIMEOUT_SECONDS)
        output = process.stdout + process.stderr
        returncode = process.returncode
        timed_out = False
    except subprocess.TimeoutExpired as error:
        output = (error.stdout or "") + (error.stderr or "")
        returncode = None
        timed_out = True
    elapsed = time.monotonic() - started
    log_path.write_text(output)
    tensor_match = re.search(r"^TENSOR_3_BY_27=(.+)$", output, re.MULTILINE)
    tensor = json.loads(tensor_match.group(1)) if tensor_match else None
    source_path = source_tensor_path(discriminant)
    source_tensor = json.loads(source_path.read_text())["tensor_3_by_27"]
    external_match = tensor is not None and tensor == source_tensor
    entries = output.count("AC1=PASS")
    verified = (
        not timed_out and returncode == 0 and entries == 18
        and "SHUFFLE_IDENTITIES=PASS" in output
        and "EXPECTED_TENSOR_MATCH=PASS" in output
        and "CERTIFICATE VERIFIED" in output
        and tensor is not None and external_match)
    return {
        "discriminant": discriminant,
        "status": "VERIFIED" if verified else "UNDECIDED",
        "seconds": elapsed,
        "entries_verified": entries,
        "entries_expected": 18,
        "shuffle_identities": "PASS" if "SHUFFLE_IDENTITIES=PASS" in output else "UNDECIDED",
        "embedded_expected_tensor_match": (
            "PASS" if "EXPECTED_TENSOR_MATCH=PASS" in output else "UNDECIDED"),
        "source_tensor_match": "PASS" if external_match else "UNDECIDED",
        "source_tensor": str(source_path.relative_to(ROOT)),
        "tensor_3_by_27": tensor,
        "returncode": returncode,
        "timed_out": timed_out,
        "run": run,
        "log": str(log_path.relative_to(ROOT)),
        "certificate": str(certificate.relative_to(ROOT)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=int, action="append")
    parser.add_argument(
        "--keep-going", action="store_true",
        help="record an UNDECIDED field and continue instead of stopping")
    args = parser.parse_args()
    wanted = {abs(value) for value in args.field} if args.field else None
    directories = sorted(CERTIFICATES.glob("*/K-*-p3"))
    directories = [directory for directory in directories
                   if (directory / "certificate.gp").exists()
                   and (wanted is None
                        or int(directory.name.split("-")[1]) in wanted)]
    if wanted is not None and len(directories) != len(wanted):
        found = {int(directory.name.split("-")[1]) for directory in directories}
        parser.error(f"missing certificate(s): {sorted(wanted - found)}")

    records = []
    for directory in directories:
        record = verify(directory)
        records.append(record)
        print(f"D={record['discriminant']:<11} {record['status']:<9} "
              f"entries={record['entries_verified']}/18 "
              f"tensor={record['source_tensor_match']} "
              f"seconds={record['seconds']:.3f}", flush=True)
        # A mismatch is a mathematical finding under this protocol.  Do not
        # continue and obscure it with later formatting or finite checks.
        # With --keep-going the finding is still recorded and the sweep
        # continues, so that one slow field does not hide the list of
        # remaining ones.
        if record["status"] != "VERIFIED" and not args.keep_going:
            break

    result_path = next_result()
    result_path.write_text(json.dumps({
        "certificates_requested": len(directories),
        "certificates_verified": sum(row["status"] == "VERIFIED"
                                     for row in records),
        "all_verified": (len(records) == len(directories)
                         and all(row["status"] == "VERIFIED" for row in records)),
        "records": records,
    }, indent=2) + "\n")
    print(f"RESULT={result_path.relative_to(ROOT)}")
    return 0 if len(records) == len(directories) and all(
        row["status"] == "VERIFIED" for row in records) else 1


if __name__ == "__main__":
    sys.exit(main())
