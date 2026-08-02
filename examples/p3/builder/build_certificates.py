#!/usr/bin/env python3
"""Run the deterministic p=3 oracle search and emit format-2 certificates.

The relative BNF is used only to find witnesses.  ``build_certificate.gp``
performs the exact arithmetic audits before it writes an entry.  Existing log
files are never opened for writing; repeated runs receive a new run number.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import select
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
TREE = HERE.parent
ROOT = HERE.parents[2]
GP_BUILDER = HERE / "build_certificate.gp"


@dataclass(frozen=True)
class Field:
    discriminant: int
    class_group: tuple[int, ...]
    tensor_path: Path


def field(discriminant: int, class_group: tuple[int, ...], path: str) -> Field:
    del path
    return Field(discriminant, class_group,
                 TREE / "source-tensors" / f"D-{abs(discriminant)}" / "tensor.json")


FIELDS = (
    field(-3640387, (18, 3, 3),
          "experiments/p3-arithmetic-pilot/results/tensor.json"),
    field(-4447704, (24, 6, 6),
          "experiments/p3-delta-sweep/results/4447704/analysis.json"),
    field(-53209523, (81, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/53209523/analysis.json"),
    field(-101375499, (72, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/101375499/analysis.json"),
    field(-134034647, (90, 18, 6),
          "experiments/p3-phase1-probes/results/rank1/134034647/analysis.json"),
    field(-138230347, (90, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/138230347/analysis.json"),
    field(-139272611, (207, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/139272611/analysis.json"),
    field(-147994487, (432, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/147994487/analysis.json"),
    field(-163004039, (252, 18, 3),
          "experiments/p3-phase1-probes/results/rank1/163004039/analysis.json"),
    field(-166596251, (279, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/166596251/analysis.json"),
    field(-198040904, (396, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/198040904/analysis.json"),
    field(-228404408, (180, 9, 3),
          "experiments/p3-phase1-probes/results/rank1/228404408/analysis.json"),
)


def expected_tensor(record: Field) -> list[list[int]]:
    data = json.loads(record.tensor_path.read_text())
    tensor = data.get("tensor_3_by_27")
    if tensor is None:
        # The legacy delta-sweep analysis used the tensor but omitted it from
        # analysis.json.  Recreate the exact object from that run's persisted
        # matrices.tsv and preserve a deterministic tensor.json snapshot.
        matrix_path = record.tensor_path.with_name("matrices.tsv")
        with matrix_path.open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        primary: dict[str, list[list[int]]] = {}
        for row in rows:
            if row["doubled"] == "1":
                continue
            matrix = primary.setdefault(row["label"], [[0] * 3 for _ in range(3)])
            column = int(row["input"]) - 1
            for output in range(3):
                matrix[output][column] = int(row[f"d{output + 1}"])

        def combine(*terms: tuple[int, list[list[int]]]) -> list[list[int]]:
            return [[sum(coefficient * matrix[i][j]
                         for coefficient, matrix in terms) % 3
                     for j in range(3)] for i in range(3)]

        d1, d2, d3 = primary["x1"], primary["x2"], primary["x3"]
        b12 = combine((1, d1), (1, d2), (-1, primary["x1+x2"]))
        b13 = combine((1, d1), (1, d3), (-1, primary["x1+x3"]))
        b23 = combine((1, d1), (1, d2), (1, d3), (-1, b12), (-1, b13),
                      (-1, primary["x1+x2+x3"]))
        tensor = [[0] * 27 for _ in range(3)]
        word = lambda i, j, k: (i * 3 + j) * 3 + k
        for i, diagonal in enumerate((d1, d2, d3)):
            for middle in range(3):
                for relation in range(3):
                    tensor[relation][word(i, middle, i)] = diagonal[middle][relation]
        for (i, k), contraction in (((0, 1), b12), ((0, 2), b13),
                                    ((1, 2), b23)):
            for middle in range(3):
                for relation in range(3):
                    value = contraction[middle][relation]
                    tensor[relation][word(i, middle, k)] = value
                    tensor[relation][word(k, middle, i)] = value
        snapshot = TREE / "source-tensors" / f"D-{abs(record.discriminant)}" / "tensor.json"
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps({
            "discriminant": record.discriminant,
            "tensor_3_by_27": tensor,
            "provenance": str(matrix_path.relative_to(ROOT)),
            "note": "Deterministic materialization of the tensor used by the legacy analysis, which omitted the tensor array from analysis.json.",
        }, indent=2) + "\n"
        if snapshot.exists() and snapshot.read_text() != payload:
            raise RuntimeError(f"legacy tensor snapshot changed: {snapshot}")
        if not snapshot.exists():
            snapshot.write_text(payload)
    if (len(tensor) != 3 or any(len(row) != 27 for row in tensor)
            or any(value not in (0, 1, 2) for row in tensor for value in row)):
        raise ValueError(f"{record.tensor_path}: tensor is not 3 by 27 over F_3")
    if data.get("discriminant", record.discriminant) != record.discriminant:
        raise ValueError(f"{record.tensor_path}: discriminant mismatch")
    return tensor


def next_run(field_dir: Path) -> tuple[int, Path, Path]:
    number = 1
    while True:
        run_dir = field_dir / "build-data" / f"run-{number:03d}"
        log_path = field_dir / f"build-{number:03d}.log"
        if not run_dir.exists() and not log_path.exists():
            run_dir.mkdir(parents=True)
            return number, run_dir, log_path
        number += 1


def stream_process(command: list[str], env: dict[str, str], log_path: Path) -> int:
    with log_path.open("x", encoding="utf-8") as log:
        process = subprocess.Popen(
            command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        last_output = time.monotonic()
        while process.poll() is None:
            ready, _, _ = select.select([process.stdout], [], [], 15.0)
            if ready:
                line = process.stdout.readline()
                if line:
                    print(line, end="", flush=True)
                    log.write(line)
                    log.flush()
                    last_output = time.monotonic()
            elif time.monotonic() - last_output >= 15.0:
                heartbeat = "BUILDER_HEARTBEAT elapsed_without_output_s=15\n"
                print(heartbeat, end="", flush=True)
                log.write(heartbeat)
                log.flush()
                last_output = time.monotonic()
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
        return process.wait()


def install_deterministic(candidate: Path, destination: Path) -> str:
    payload = candidate.read_bytes()
    if destination.exists():
        if destination.read_bytes() != payload:
            raise RuntimeError(
                f"determinism failure: {destination} differs from rebuilt certificate")
        candidate.unlink()
        return "UNCHANGED"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    candidate.unlink()
    return "CREATED"


def build_one(record: Field) -> dict[str, object]:
    name = f"K-{abs(record.discriminant)}-p3"
    field_dir = TREE / "certificates" / name
    run, run_dir, log_path = next_run(field_dir)
    candidate = run_dir / "certificate.gp.tmp"
    tensor = expected_tensor(record)
    environment = os.environ.copy()
    environment.update({
        "P3_DISC": str(record.discriminant),
        "P3_EXPECTED_CYC": str(list(record.class_group)),
        "P3_RESULT_DIR": str(run_dir),
        "P3_CERT_PATH": str(candidate),
        "P3_EXPECTED_TENSOR": json.dumps(tensor, separators=(",", ":")),
    })
    started = time.monotonic()
    print(f"BUILD {name} run={run:03d}", flush=True)
    returncode = stream_process(["gp", "-qf", str(GP_BUILDER)],
                                environment, log_path)
    elapsed = time.monotonic() - started
    log_text = log_path.read_text()
    completed = "ARITHMETIC PILOT COMPLETE" in log_text
    gp_error = any(
        "***" in line and "Warning:" not in line
        for line in log_text.splitlines())
    if returncode or not candidate.exists() or not completed or gp_error:
        return {
            "discriminant": record.discriminant,
            "status": "UNDECIDED",
            "reason": (
                f"builder exit status {returncode}; "
                f"completion_marker={completed}; gp_error={gp_error}"),
            "seconds": elapsed,
            "run": run,
            "log": str(log_path.relative_to(ROOT)),
        }
    status = install_deterministic(candidate, field_dir / "certificate.gp")
    return {
        "discriminant": record.discriminant,
        "status": status,
        "seconds": elapsed,
        "run": run,
        "log": str(log_path.relative_to(ROOT)),
        "certificate": str((field_dir / "certificate.gp").relative_to(ROOT)),
        "source_tensor": str(record.tensor_path.relative_to(ROOT)),
    }


def unique_results_path(stem: str) -> Path:
    results = TREE / "results"
    results.mkdir(exist_ok=True)
    number = 1
    while (results / f"{stem}-{number:03d}.json").exists():
        number += 1
    return results / f"{stem}-{number:03d}.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--field", type=int, action="append",
        help="absolute or signed discriminant; default: all twelve")
    args = parser.parse_args()
    wanted = {abs(value) for value in args.field} if args.field else None
    selected = [record for record in FIELDS
                if wanted is None or abs(record.discriminant) in wanted]
    if wanted is not None and len(selected) != len(wanted):
        known = {abs(record.discriminant) for record in selected}
        parser.error(f"unknown field(s): {sorted(wanted - known)}")

    records = []
    for record in selected:
        result = build_one(record)
        records.append(result)
        print(json.dumps(result, sort_keys=True), flush=True)
        if result["status"] == "UNDECIDED":
            break

    output = unique_results_path("build-costs")
    output.write_text(json.dumps({
        "randomness_used": False,
        "seeds": [],
        "fields_requested": len(selected),
        "fields_completed": len(records),
        "records": records,
    }, indent=2) + "\n")
    return 1 if any(row["status"] == "UNDECIDED" for row in records) else 0


if __name__ == "__main__":
    sys.exit(main())
