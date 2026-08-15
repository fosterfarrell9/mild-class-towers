#!/usr/bin/env python3
"""Recheck strong-freeness certificates on C-verifier tensors."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(HERE)]

from admissible import find_anick_witness  # noqa: E402
from extension_search import FiniteField, coset_anick_search  # noqa: E402
from hilbert_depth import measured_comparison  # noqa: E402

PILOT_DISC = -3640387
GROEBNER_DISC = -4447704
F9_DISC = -139272611
RATIONAL_ANICK = {
    -53209523, -101375499, -134034647, -138230347, -147994487,
    -163004039, -166596251, -198040904, -228404408,
}


def latest_verification() -> Path:
    """Newest verification of the twelve fields these checks are about.

    The results directory also holds the block verification of all 2497
    fields; picking the newest file outright would find that one and fail
    the count test below.
    """
    paths = sorted((HERE.parent / "results").glob("verification-*.json"))
    for path in reversed(paths):
        data = json.loads(path.read_text())
        if data.get("all_verified") and data.get("certificates_verified") == 12:
            return path
    raise FileNotFoundError(
        "no complete 12/12 verification in results/; run verifier/verify_all.py")


def next_result() -> Path:
    number = 1
    results = HERE.parent / "results"
    while (results / f"strong-freeness-{number:03d}.json").exists():
        number += 1
    return results / f"strong-freeness-{number:03d}.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path)
    args = parser.parse_args()
    verification_path = (args.verification or latest_verification()).resolve()
    verification = json.loads(verification_path.read_text())
    if not verification["all_verified"] or verification["certificates_verified"] != 12:
        raise RuntimeError("strong-freeness checks require a complete 12/12 C verification")
    tensors = {record["discriminant"]: record["tensor_3_by_27"]
               for record in verification["records"]}
    expected = RATIONAL_ANICK | {PILOT_DISC, GROEBNER_DISC, F9_DISC}
    if set(tensors) != expected:
        raise RuntimeError("verification result does not contain the fixed twelve fields")

    records: list[dict[str, object]] = []
    f9 = FiniteField((1, 0, 1))
    for discriminant in sorted(tensors):
        tensor = tensors[discriminant]
        started = time.monotonic()
        if discriminant in RATIONAL_ANICK:
            witness = find_anick_witness(tensor)
            status = "STRONGLY_FREE" if witness["found"] else "UNDECIDED"
            record = {
                "discriminant": discriminant,
                "certificate": "ANICK_F3",
                "status": status,
                "witness": witness,
            }
        elif discriminant == F9_DISC:
            witness = coset_anick_search(tensor, f9)
            valid = (witness["found"] and witness["exhaustive"]
                     and witness["tested"] == 663390)
            record = {
                "discriminant": discriminant,
                "certificate": "ANICK_F9_WITH_F3_DESCENT",
                "status": "STRONGLY_FREE" if valid else "UNDECIDED",
                "witness": witness,
            }
        elif discriminant == GROEBNER_DISC:
            witness = measured_comparison(tensor, 10, f"D={discriminant}")
            valid = (witness["verdict"] == "STRONGLY_FREE"
                     and witness["terminated_completion"]
                     and witness["positive_certificate_comparison_bound"] == 25)
            record = {
                "discriminant": discriminant,
                "certificate": "TERMINATING_GROEBNER_AUTOMATON",
                "status": "STRONGLY_FREE" if valid else "UNDECIDED",
                "witness": witness,
            }
        else:
            # The pilot has no positive certificate.  Exhaust both finite
            # Anick searches again; their failure is evidence only and is not
            # promoted to a negative strong-freeness statement.
            rational = find_anick_witness(tensor)
            extension = coset_anick_search(tensor, f9)
            if rational["found"] or extension["found"]:
                status = "STRONGLY_FREE"
            else:
                status = "UNDECIDED"
            record = {
                "discriminant": discriminant,
                "certificate": "NO_POSITIVE_CERTIFICATE",
                "status": status,
                "rational_anick": rational,
                "f9_anick": extension,
                "interpretation": (
                    "Failure of the two Anick searches is not a negative "
                    "strong-freeness certificate."),
            }
        record["seconds"] = time.monotonic() - started
        records.append(record)
        print(f"D={discriminant:<11} {record['certificate']:<31} "
              f"{record['status']} seconds={record['seconds']:.3f}", flush=True)
        if discriminant != PILOT_DISC and record["status"] != "STRONGLY_FREE":
            break

    theorem_fields = [row for row in records if row["discriminant"] != PILOT_DISC]
    theorem_verified = sum(row["status"] == "STRONGLY_FREE"
                           for row in theorem_fields)
    complete = len(theorem_fields) == 11 and theorem_verified == 11
    output = next_result()
    output.write_text(json.dumps({
        "source_verification": str(verification_path.relative_to(ROOT)),
        "randomness_used": False,
        "seeds": [],
        "theorem_fields_expected": 11,
        "theorem_fields_strongly_free": theorem_verified,
        "theorem_chain_complete": complete,
        "pilot_status": next(row["status"] for row in records
                             if row["discriminant"] == PILOT_DISC),
        "records": records,
    }, indent=2) + "\n")
    print(f"RESULT={output.relative_to(ROOT)}")
    return 0 if complete else 1


if __name__ == "__main__":
    sys.exit(main())
