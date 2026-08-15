#!/usr/bin/env python3
"""Negative tests for the mandatory verifier rejection paths.

The unified verifier serves every odd prime of the repository, so the
suite mutates one fixture per certificate dialect: the p = 3 fixture
(x-labels, embedded expected tensor, J-corrected classes) and the
p = 5 fixture (letter labels, no embedded tensor).  For each fixture
the five certificate mutations of the paper's verification appendix
are replayed --- the character vector, the normalized automorphism,
the multiplicity of an entry, a norm-class vector, and the absolute
field model of a temporary copy are altered separately --- and the
verifier must reject every copy with the expected message; the
remaining cases cover structural rejections of the certificate
container.  `make check` runs the unmodified fixtures of all three
primes first and then this script; the outcomes are recorded in
`rejection-tests.json`.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFIER = HERE / "verify_certificate"
SOURCE_P3 = (HERE.parent / "certificates" / "p3"
             / "000" / "K-3640387-p3" / "certificate.gp")
SOURCE_P5 = HERE.parent / "certificates" / "p5" / "K-2800905-p5" / "certificate.gp"

# The only length-3 integral ~-vectors of an entry line are the
# character vector and the norm-class vector, in this order; the
# stored sigma of the p = 5 fixture is the first length-10 integral
# ~-vector of an entry.
THREE_VECTOR = re.compile(r"\[(-?\d+), (-?\d+), (-?\d+)\]~")
TEN_VECTOR = re.compile(r"\[-?\d+(?:, -?\d+){9}\]~")


def run(path: Path) -> tuple[int, str]:
    process = subprocess.run([str(VERIFIER), str(path)], capture_output=True,
                             text=True, timeout=300)
    return process.returncode, process.stdout + process.stderr


def p3_cases() -> dict[str, tuple[str, str]]:
    source = SOURCE_P3.read_text()
    lines = source.splitlines(keepends=True)
    if len(lines) != 20 or not lines[-2].startswith('["x1+x3", 3,'):
        raise RuntimeError("unexpected p3 fixture line layout")
    if not (lines[1].startswith('["x1", 1,')
            and lines[2].startswith('["x1", 2,')):
        raise RuntimeError("unexpected p3 fixture entry order")

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
    cases["format"] = (
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

    return cases


def p5_cases() -> dict[str, tuple[str, str]]:
    source = SOURCE_P5.read_text()
    lines = source.splitlines(keepends=True)
    if len(lines) != 18 or not lines[-1].startswith('["b+c",3,'):
        raise RuntimeError("unexpected p5 fixture line layout")
    if not lines[1].startswith('["a",2,'):
        raise RuntimeError("unexpected p5 fixture entry order")

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
        "".join(shortened), "certificate must contain exactly 18 entries")

    return cases


def main() -> int:
    suites = {"p3": p3_cases(), "p5": p5_cases()}
    records = []
    with tempfile.TemporaryDirectory(prefix="cert-negative-") as temp:
        temp_dir = Path(temp)
        for fixture, cases in suites.items():
            for name, (payload, expected_error) in cases.items():
                path = temp_dir / f"{fixture}-{name}.gp"
                path.write_text(payload)
                returncode, output = run(path)
                passed = returncode != 0 and expected_error in output
                records.append({
                    "fixture": fixture,
                    "case": name,
                    "status": "PASS" if passed else "FAIL",
                    "expected_error": expected_error,
                    "returncode": returncode,
                })
                print(f"{fixture}/{name}: {'PASS' if passed else 'FAIL'}")

    result = HERE / "rejection-tests.json"
    result.write_text(json.dumps({
        "all_passed": all(row["status"] == "PASS" for row in records),
        "records": records,
    }, indent=2) + "\n")
    return 0 if all(row["status"] == "PASS" for row in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
