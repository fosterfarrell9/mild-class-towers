#!/usr/bin/env python3
"""Small negative tests for mandatory p=3 verifier rejection paths."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFIER = HERE / "verify_certificate"
SOURCE = HERE.parent / "certificates" / "K-3640387-p3" / "certificate.gp"


def run(path: Path) -> tuple[int, str]:
    process = subprocess.run([str(VERIFIER), str(path)], capture_output=True,
                             text=True, timeout=30)
    return process.returncode, process.stdout + process.stderr


def main() -> int:
    source = SOURCE.read_text()
    cases: dict[str, tuple[str, str]] = {}
    cases["format_1"] = (
        source.replace("[2,", "[1,", 1), "unsupported certificate format")

    lines = source.splitlines(keepends=True)
    if len(lines) != 20 or not lines[-2].startswith('["x1+x3", 3,'):
        raise RuntimeError("unexpected fixture line layout")
    lines[-3] = lines[-3].rstrip("\n")
    if not lines[-3].endswith(","):
        raise RuntimeError("penultimate entry has no separator")
    lines[-3] = lines[-3][:-1] + "\n"
    del lines[-2]
    cases["missing_entry"] = (
        "".join(lines), "must contain exactly 18 entries")

    needle = "[[0, 0, 0, 0, 2, 2, 0, 1, 0"
    if source.count(needle) != 1:
        raise RuntimeError("expected-tensor fixture marker is not unique")
    cases["wrong_expected_tensor"] = (
        source.replace(needle, "[[1, 0, 0, 0, 2, 2, 0, 1, 0", 1),
        "reconstructed tensor disagrees with expected tensor")

    records = []
    with tempfile.TemporaryDirectory(prefix="p3-cert-negative-") as temp:
        temp_dir = Path(temp)
        for name, (payload, expected_error) in cases.items():
            path = temp_dir / f"{name}.gp"
            path.write_text(payload)
            returncode, output = run(path)
            passed = returncode != 0 and expected_error in output
            records.append({
                "case": name,
                "status": "PASS" if passed else "FAIL",
                "expected_error": expected_error,
                "returncode": returncode,
            })
            print(f"{name}: {'PASS' if passed else 'FAIL'}")

    result = HERE.parent / "results" / "rejection-tests.json"
    result.write_text(json.dumps({
        "all_passed": all(row["status"] == "PASS" for row in records),
        "records": records,
    }, indent=2) + "\n")
    return 0 if all(row["status"] == "PASS" for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
