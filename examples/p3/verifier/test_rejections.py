#!/usr/bin/env python3
"""Negative tests for the mandatory verifier rejection paths.

The first five cases replay the certificate mutations of the paper's
verification appendix: the character vector, the normalized
automorphism, the multiplicity of an entry, a norm-class vector, and
the absolute field model of a temporary copy are altered separately,
and the verifier must reject every copy with the expected message.
The remaining cases cover structural rejections of the certificate
container.  `make check` runs the unmodified fixture first and then
this script; the outcomes are recorded in
`../results/rejection-tests.json`.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFIER = HERE / "verify_certificate"
SOURCE = HERE.parent / "certificates" / "K-3640387-p3" / "certificate.gp"

# The only length-3 integral ~-vectors of an entry line are the
# character vector and the norm-class vector, in this order.
THREE_VECTOR = re.compile(r"\[(-?\d+), (-?\d+), (-?\d+)\]~")


def run(path: Path) -> tuple[int, str]:
    process = subprocess.run([str(VERIFIER), str(path)], capture_output=True,
                             text=True, timeout=120)
    return process.returncode, process.stdout + process.stderr


def main() -> int:
    source = SOURCE.read_text()
    lines = source.splitlines(keepends=True)
    if len(lines) != 20 or not lines[-2].startswith('["x1+x3", 3,'):
        raise RuntimeError("unexpected fixture line layout")
    if not (lines[1].startswith('["x1", 1,')
            and lines[2].startswith('["x1", 2,')):
        raise RuntimeError("unexpected fixture entry order")

    cases: dict[str, tuple[str, str]] = {}

    # 1. The character vector: [1, 0, 0]~ of entry x1/e1 no longer
    #    matches the label x1.
    needle = '["x1", 1, [1, 0, 0]~,'
    if source.count(needle) != 1:
        raise RuntimeError("character-vector fixture marker is not unique")
    cases["character_vector"] = (
        source.replace(needle, '["x1", 1, [1, 1, 0]~,', 1),
        "character label/vector mismatch")

    # 2. The normalized automorphism: one coefficient of the stored
    #    sigma of entry x1/e1 is changed.
    sigma = "[21283, -42566, -21283, 35162, 35162, 0]~"
    if source.count(sigma) != 3:
        raise RuntimeError("sigma fixture marker is not shared by the x1 entries")
    cases["normalized_automorphism"] = (
        source.replace(sigma, "[21284, -42566, -21283, 35162, 35162, 0]~", 1),
        "stored sigma does not fix K")

    # 3. The multiplicity of an entry: the line of x1/e2 is replaced
    #    by a second copy of x1/e1, so the count stays at 18.
    duplicated = lines.copy()
    duplicated[2] = duplicated[1]
    cases["entry_multiplicity"] = (
        "".join(duplicated), "duplicate certificate entry")

    # 4. A norm-class vector: the second length-3 vector of entry
    #    x1/e1 is its stored norm class; one coordinate is moved.
    matches = list(THREE_VECTOR.finditer(lines[1]))
    if len(matches) != 2:
        raise RuntimeError("entry x1/e1 does not carry exactly two 3-vectors")
    last = matches[-1]
    first_coordinate = int(last.group(1))
    if not 0 <= first_coordinate <= 2:
        raise RuntimeError("stored norm class is not reduced modulo 3")
    mutated_entry = (
        lines[1][:last.start()]
        + "[%d, %s, %s]~" % ((first_coordinate + 1) % 3,
                             last.group(2), last.group(3))
        + lines[1][last.end():])
    mutated = lines.copy()
    mutated[1] = mutated_entry
    cases["norm_class_vector"] = (
        "".join(mutated), "J-corrected class coordinates mismatch")

    # 5. The absolute field model: the constant term of the stored
    #    absolute polynomial of the x1 entries is changed in x1/e1.
    f_abs = ("x^6 + 96*x^4 + 2304*x^2 + "
             "2648356114460005728290905046600437507")
    if source.count(f_abs) != 3:
        raise RuntimeError("absolute-model fixture marker is not shared")
    cases["absolute_field_model"] = (
        source.replace(f_abs, f_abs[:-1] + "9", 1),
        "relative/absolute field models are incompatible")

    # Structural rejections of the certificate container.
    cases["format_1"] = (
        source.replace("[2,", "[1,", 1), "unsupported certificate format")

    shortened = lines.copy()
    shortened[-3] = shortened[-3].rstrip("\n")
    if not shortened[-3].endswith(","):
        raise RuntimeError("penultimate entry has no separator")
    shortened[-3] = shortened[-3][:-1] + "\n"
    del shortened[-2]
    cases["missing_entry"] = (
        "".join(shortened), "must contain exactly 18 entries")

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
