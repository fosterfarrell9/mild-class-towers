#!/usr/bin/env python3
"""Negative tests for the mandatory p=5 verifier rejection paths.

The first five cases replay the certificate mutations of the paper's
verification appendix: the character vector, the normalized
automorphism, the multiplicity of an entry, a norm-class vector, and
the absolute field model of a temporary copy are altered separately,
and the verifier must reject every copy with the expected message.
The remaining cases cover structural rejections of the certificate
container; the tensor comparison has no stored-tensor case here
because the p=5 certificate carries no expected tensor --- that
comparison runs against a separate result record.  `make check` runs
the unmodified fixture first and then this script; the outcomes are
recorded in `rejection-tests.json`.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFIER = HERE / "verify_certificate"
SOURCE = HERE / "K-2800905-p5" / "certificate.gp"

# The only length-3 integral ~-vectors of an entry are the character
# vector and the norm-class vector, in this order; the stored sigma
# is the first length-10 integral ~-vector of an entry.
THREE_VECTOR = re.compile(r"\[(-?\d+), (-?\d+), (-?\d+)\]~")
TEN_VECTOR = re.compile(r"\[-?\d+(?:, -?\d+){9}\]~")


def run(path: Path) -> tuple[int, str]:
    process = subprocess.run([str(VERIFIER), str(path)], capture_output=True,
                             text=True, timeout=120)
    return process.returncode, process.stdout + process.stderr


def main() -> int:
    source = SOURCE.read_text()
    lines = source.splitlines(keepends=True)
    if len(lines) != 18 or not lines[-1].startswith('["b+c",3,'):
        raise RuntimeError("unexpected fixture line layout")
    if not lines[1].startswith('["a",2,'):
        raise RuntimeError("unexpected fixture entry order")

    cases: dict[str, tuple[str, str]] = {}

    # 1. The character vector: [1, 0, 0]~ of entry a/1 no longer
    #    matches the label a.
    needle = '["a",1,[1, 0, 0]~,'
    if source.count(needle) != 1:
        raise RuntimeError("character-vector fixture marker is not unique")
    cases["character_vector"] = (
        source.replace(needle, '["a",1,[1, 1, 0]~,', 1),
        "character label/vector mismatch")

    # 2. The normalized automorphism: the first coefficient of the
    #    stored sigma of entry a/1 (the first length-10 vector of the
    #    first line) is changed.
    sigma = TEN_VECTOR.search(lines[0])
    if sigma is None:
        raise RuntimeError("entry a/1 carries no stored sigma")
    bumped = re.sub(r"^\[(-?\d+)",
                    lambda m: "[%d" % (int(m.group(1)) + 1),
                    sigma.group(0), count=1)
    mutated = lines.copy()
    mutated[0] = lines[0][:sigma.start()] + bumped + lines[0][sigma.end():]
    cases["normalized_automorphism"] = (
        "".join(mutated), "stored sigma does not fix K")

    # 3. The multiplicity of an entry: the entry a/2 is relabelled as
    #    a second a/1, so the count stays at 18.
    needle = '["a",2,'
    if source.count(needle) != 1:
        raise RuntimeError("multiplicity fixture marker is not unique")
    cases["entry_multiplicity"] = (
        source.replace(needle, '["a",1,', 1), "duplicate certificate entry")

    # 4. A norm-class vector: the second length-3 vector of the first
    #    line is the stored norm class of entry a/1; one coordinate
    #    is moved.
    matches = list(THREE_VECTOR.finditer(lines[0]))
    if len(matches) != 2:
        raise RuntimeError("entry a/1 does not carry exactly two 3-vectors")
    last = matches[-1]
    first_coordinate = int(last.group(1))
    if not 0 <= first_coordinate <= 4:
        raise RuntimeError("stored norm class is not reduced modulo 5")
    mutated = lines.copy()
    mutated[0] = (
        lines[0][:last.start()]
        + "[%d, %s, %s]~" % ((first_coordinate + 1) % 5,
                             last.group(2), last.group(3))
        + lines[0][last.end():])
    cases["norm_class_vector"] = (
        "".join(mutated), "norm-class coordinates mismatch")

    # 5. The absolute field model: the constant term of the stored
    #    absolute polynomial of the a entries is changed in a/1.
    f_abs = ("y^10 - 2*y^9 - 2053*y^8 - 65684*y^7 - 174901*y^6 + "
             "64171174*y^5 + 2665457137*y^4 + 66051314232*y^3 + "
             "1091726153376*y^2 + 9461326049280*y + 31436965969920")
    if source.count(f_abs) != 3:
        raise RuntimeError("absolute-model fixture marker is not shared")
    cases["absolute_field_model"] = (
        source.replace(f_abs, f_abs[:-1] + "2", 1),
        "relative/absolute field models are incompatible")

    # Structural rejections of the certificate container.
    cases["format"] = (
        source.replace("[2,", "[1,", 1), "unsupported certificate format")

    shortened = lines.copy()
    if not shortened[1].rstrip().endswith("],"):
        raise RuntimeError("entry a/2 does not close with a separator")
    del shortened[1]
    cases["missing_entry"] = (
        "".join(shortened), "certificate must contain 18 or 27 entries")

    records = []
    with tempfile.TemporaryDirectory(prefix="p5-cert-negative-") as temp:
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

    result = HERE / "rejection-tests.json"
    result.write_text(json.dumps({
        "all_passed": all(row["status"] == "PASS" for row in records),
        "records": records,
    }, indent=2) + "\n")
    return 0 if all(row["status"] == "PASS" for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
