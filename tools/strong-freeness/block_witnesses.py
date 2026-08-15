#!/usr/bin/env python3
"""Strong-freeness witnesses for the block fields beyond the criterion.

The cone criterion covers eighty-seven fields of the block from the
verified tensors alone.  Five further fields carry a direct
strong-freeness witness: two have a rational Anick witness at a
degenerate cone point, and three have a terminating Gröbner
completion.  This driver recomputes all five witnesses from the
C-verifier tensors of the latest complete block verification and
writes ``results/block-witnesses-<n>.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE)]

from admissible import find_anick_witness  # noqa: E402
from hilbert_depth import measured_comparison  # noqa: E402

RATIONAL_ANICK = (-211248887, -263780072)
GROEBNER = (-4447704, -192928619, -263310215)
BLOCK_SIZE = 2497


def latest_verification() -> Path:
    """Newest complete block verification.

    Selected by size, not by date: the results directory also holds the
    twelve-field verifications of the theorem fields.
    """
    paths = sorted((ROOT / "records" / "p3" / "results").glob("verification-*.json"))
    for path in reversed(paths):
        data = json.loads(path.read_text())
        if (data.get("all_verified")
                and data.get("certificates_verified") == BLOCK_SIZE):
            return path
    raise FileNotFoundError(
        "no complete block verification in records/p3/results/; run tools/verify_all_p3.py")


def next_result() -> Path:
    number = 1
    results = ROOT / "records" / "p3" / "results"
    while (results / f"block-witnesses-{number:03d}.json").exists():
        number += 1
    return results / f"block-witnesses-{number:03d}.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path)
    args = parser.parse_args()
    verification_path = (args.verification or latest_verification()).resolve()
    verification = json.loads(verification_path.read_text())
    if (not verification["all_verified"]
            or verification["certificates_verified"] != BLOCK_SIZE):
        raise RuntimeError(
            "block witnesses require a complete block C verification")
    tensors = {record["discriminant"]: record["tensor_3_by_27"]
               for record in verification["records"]}

    records: list[dict[str, object]] = []
    for discriminant in RATIONAL_ANICK:
        started = time.monotonic()
        witness = find_anick_witness(tensors[discriminant])
        records.append({
            "discriminant": discriminant,
            "certificate": "ANICK_F3",
            "status": "STRONGLY_FREE" if witness["found"] else "UNDECIDED",
            "witness": witness,
            "seconds": time.monotonic() - started,
        })
    for discriminant in GROEBNER:
        started = time.monotonic()
        witness = measured_comparison(tensors[discriminant], 10,
                                      f"D={discriminant}")
        valid = (witness["verdict"] == "STRONGLY_FREE"
                 and witness["terminated_completion"]
                 and witness["positive_certificate_comparison_bound"]
                 is not None)
        records.append({
            "discriminant": discriminant,
            "certificate": "TERMINATING_GROEBNER",
            "status": "STRONGLY_FREE" if valid else "UNDECIDED",
            "witness": witness,
            "seconds": time.monotonic() - started,
        })
    for record in records:
        print(f"D={record['discriminant']:<11} {record['certificate']:<22} "
              f"{record['status']} seconds={record['seconds']:.3f}",
              flush=True)

    output = next_result()
    output.write_text(json.dumps({
        "verification": verification_path.name,
        "randomness_used": False,
        "seeds": [],
        "strongly_free": sum(r["status"] == "STRONGLY_FREE"
                             for r in records),
        "of": len(records),
        "records": records,
    }, indent=2) + "\n")
    print(f"RESULT={output.relative_to(ROOT)}")
    return 0 if all(r["status"] == "STRONGLY_FREE" for r in records) else 1


if __name__ == "__main__":
    sys.exit(main())
