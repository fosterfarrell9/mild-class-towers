#!/usr/bin/env python3
"""Exhaustive rational Anick search over the criterion-undecided fields.

For every p=3 field on which the cone criterion is silent, this driver
runs ``admissible.find_anick_witness`` on the verified tensor: an
exhaustive search over all 11232 elements of GL_3(F_3) for a rational
change of variables whose three cubic high terms are combinatorially
free.  A found witness proves the cubic relation space strongly free;
an exhausted search proves that no rational Anick witness exists.

The tensors are read from the complete verification records under
``records/p3/results/``, the criterion verdicts from
``records/p3/cone-criterion/``.  Nothing is recomputed and nothing in
the repository is written; every shard writes one JSON file to the
directory named by ``--out``.

Sharded run (24 processes), then merge:

    python3 tools/strong-freeness/anick_scan.py --shards 24 --shard 0 \
        --out ~/anick-scan/shard-00.json
    ...
    python3 tools/strong-freeness/anick_scan.py \
        --merge '~/anick-scan/shard-*.json' \
        --out ~/anick-scan/anick-scan-merged.json
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE)]

from admissible import find_anick_witness  # noqa: E402


def tensor_rank_mod3(tensor: list[list[int]]) -> int:
    rows = [[value % 3 for value in row] for row in tensor]
    rank = 0
    for column in range(27):
        pivot = next((r for r in range(rank, len(rows))
                      if rows[r][column] % 3), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = 1 if rows[rank][column] % 3 == 1 else 2
        rows[rank] = [(value * inverse) % 3 for value in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][column] % 3:
                factor = rows[r][column]
                rows[r] = [(rows[r][i] - factor * rows[rank][i]) % 3
                           for i in range(27)]
        rank += 1
        if rank == len(rows):
            break
    return rank

TENSOR_SOURCES = (
    "records/p3/results/verification.json",
    "records/p3/results/verification-001.json",
    "records/p3/results/verification-002.json",
    "records/p3/results/verification-003.json",
    "records/p3/results/verification-004.json",
    "records/p3/results/verification-005.json",
    "records/p3/results/block23-verification-starthinker.json",
    "records/p3/results/slice2/verification-starthinker.json",
)

CRITERION_SOURCES = (
    "records/p3/cone-criterion/summary.json",
    "records/p3/cone-criterion/summary-block23-starthinker.json",
    "records/p3/cone-criterion/summary-slice2.json",
)

EXPECTED_TARGETS = 12244
GL3_SIZE = 11232


def load_tensors() -> dict[int, list[list[int]]]:
    tensors: dict[int, list[list[int]]] = {}
    for rel in TENSOR_SOURCES:
        data = json.loads((ROOT / rel).read_text())
        records = data.get("records", data) if isinstance(data, dict) else data
        for record in records:
            if isinstance(record, dict) and "tensor_3_by_27" in record:
                tensors[record["discriminant"]] = record["tensor_3_by_27"]
    return tensors


def load_targets() -> list[int]:
    undecided: set[int] = set()
    seen: set[int] = set()
    for rel in CRITERION_SOURCES:
        data = json.loads((ROOT / rel).read_text())
        records = data.get("records", data) if isinstance(data, dict) else data
        for record in records:
            discriminant = record["discriminant"]
            if discriminant in seen:
                raise RuntimeError(f"duplicate criterion record {discriminant}")
            seen.add(discriminant)
            if record["min_transverse_cone_degree"] is None:
                undecided.add(discriminant)
    if len(undecided) != EXPECTED_TARGETS:
        raise RuntimeError(
            f"expected {EXPECTED_TARGETS} criterion-undecided fields, "
            f"found {len(undecided)} in {len(seen)} records")
    return sorted(undecided, key=abs)


def scan(targets: list[int], tensors: dict[int, list[list[int]]],
         label: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    started = time.monotonic()
    for index, discriminant in enumerate(targets, 1):
        field_started = time.monotonic()
        rank = tensor_rank_mod3(tensors[discriminant])
        if rank < 3:
            records.append({
                "discriminant": discriminant,
                "found": False,
                "tested": 0,
                "tensor_rank": rank,
                "note": "tensor rank below three, no Anick witness possible",
                "seconds": time.monotonic() - field_started,
            })
            print(f"RANK-DEFECT D={discriminant} rank={rank}, skipped",
                  flush=True)
            continue
        witness = find_anick_witness(tensors[discriminant])
        record: dict[str, object] = {
            "discriminant": discriminant,
            "found": witness["found"],
            "tested": witness["tested"],
            "seconds": time.monotonic() - field_started,
        }
        if witness["found"]:
            record["matrix"] = witness["matrix"]
            record["leaders"] = witness["leaders"]
            print(f"WITNESS D={discriminant} after {witness['tested']} "
                  f"of {GL3_SIZE}", flush=True)
        records.append(record)
        if index % 25 == 0 or index == len(targets):
            print(f"{label}: {index}/{len(targets)} scanned, "
                  f"{sum(1 for r in records if r['found'])} witnesses, "
                  f"{time.monotonic() - started:.0f}s", flush=True)
    return records


def merge(pattern: str, out: Path) -> int:
    shards = sorted(glob.glob(str(Path(pattern).expanduser())))
    if not shards:
        raise RuntimeError(f"no shard files match {pattern}")
    by_discriminant: dict[int, dict[str, object]] = {}
    of = None
    for path in shards:
        shard = json.loads(Path(path).read_text())
        of = shard["of"]
        for record in shard["records"]:
            by_discriminant[record["discriminant"]] = record
    records = [by_discriminant[d] for d in sorted(by_discriminant, key=abs)]
    witnesses = [r for r in records if r["found"]]
    summary = {
        "engine": "python",
        "routine": "admissible.find_anick_witness",
        "gl3_size": GL3_SIZE,
        "shards": len(shards),
        "shards_expected": of,
        "fields": len(records),
        "fields_expected": EXPECTED_TARGETS,
        "witnesses_found": len(witnesses),
        "witness_discriminants": [r["discriminant"] for r in witnesses],
        "records": records,
    }
    out.expanduser().write_text(json.dumps(summary, indent=1) + "\n")
    complete = (len(records) == EXPECTED_TARGETS and len(shards) == of)
    print(f"merged {len(shards)} shards: {len(records)} fields, "
          f"{len(witnesses)} witnesses, "
          f"{'COMPLETE' if complete else 'INCOMPLETE'}", flush=True)
    return 0 if complete else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--only", type=str,
                        help="comma-separated discriminants, overrides sharding")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--merge", type=str,
                        help="glob of shard files to merge; no scanning")
    args = parser.parse_args()

    if args.merge:
        return merge(args.merge, args.out)

    tensors = load_tensors()
    targets = load_targets()
    missing = [d for d in targets if d not in tensors]
    if missing:
        raise RuntimeError(f"{len(missing)} targets without tensor, "
                           f"first {missing[:3]}")

    if args.only:
        wanted = [int(part) for part in args.only.split(",")]
        for discriminant in wanted:
            if discriminant not in set(targets):
                raise RuntimeError(f"{discriminant} is not criterion-undecided")
        selected = wanted
        label = "only"
    else:
        if not 0 <= args.shard < args.shards:
            raise RuntimeError("--shard out of range")
        selected = [d for i, d in enumerate(targets)
                    if i % args.shards == args.shard]
        label = f"shard {args.shard}/{args.shards}"

    records = scan(selected, tensors, label)
    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "engine": "python",
        "routine": "admissible.find_anick_witness",
        "gl3_size": GL3_SIZE,
        "shard": args.shard if not args.only else None,
        "of": args.shards if not args.only else None,
        "only": args.only,
        "fields": len(records),
        "records": records,
    }, indent=1) + "\n")
    print(f"RESULT={out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
