#!/usr/bin/env python3
"""Deterministic PARI cross-checks for a small MJ table validation sample."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import mj_classgroup_tables as mj


FAMILIES = ("cl3mod8", "cl7mod8", "cl4mod16", "cl8mod16")
TARGETS = (1_000_000, 6_000_000, 11_000_000, 16_000_000, 21_000_000)


def select_samples(data_dir: Path) -> list[mj.Record]:
    samples: list[mj.Record] = []
    for family in FAMILIES:
        path = data_dir / f"{family}.0.gz"
        target_index = 0
        for record in mj.read_table(path):
            magnitude = -record.discriminant
            while target_index < len(TARGETS) and magnitude >= TARGETS[target_index]:
                samples.append(record)
                target_index += 1
            if target_index == len(TARGETS):
                break
        if target_index != len(TARGETS):
            raise RuntimeError(f"{path}: did not reach all deterministic targets")
    return samples


def pari_rows(samples: list[mj.Record]) -> dict[int, tuple[int, str]]:
    commands = []
    for record in samples:
        commands.append(
            f'q=quadclassunit({record.discriminant});'
            f'print({record.discriminant},"\\t",q.no,"\\t",q.cyc)'
        )
    commands.append("quit")
    process = subprocess.run(
        ["gp", "-q"],
        input=";\n".join(commands) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        raise RuntimeError(f"GP failed:\n{process.stderr}")
    rows: dict[int, tuple[int, str]] = {}
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            raise RuntimeError(f"unexpected GP output: {line!r}")
        rows[int(fields[0])] = (int(fields[1]), fields[2].replace(" ", ""))
    return rows


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    args = parser.parse_args(argv)
    try:
        samples = select_samples(args.data_dir)
        checked = pari_rows(samples)
        lengths: set[int] = set()
        class_numbers: set[int] = set()
        for record in samples:
            pari_number, pari_cyc = checked[record.discriminant]
            expected_cyc = list(record.invariants)
            # MJ explicitly stores the trivial group as [1]; PARI uses [].
            if expected_cyc == [1]:
                expected_cyc = []
            expected_text = "[" + ",".join(str(c) for c in expected_cyc) + "]"
            if pari_number != record.class_number or pari_cyc != expected_text:
                raise RuntimeError(
                    f"D={record.discriminant}: table h={record.class_number}, "
                    f"cyc={record.invariants}; PARI h={pari_number}, cyc={pari_cyc}"
                )
            lengths.add(len(record.invariants))
            class_numbers.add(record.class_number)
            print(
                "PARI_MATCH "
                f"D={record.discriminant} "
                f"cyc={list(record.invariants)} "
                f"h={record.class_number} "
                f"source={record.source.path.name}:{record.record_number}"
            )
        if len(samples) != 20 or len(lengths) < 2 or len(class_numbers) < 2:
            raise RuntimeError("sample diversity requirement not met")
    except (OSError, mj.TableFormatError, RuntimeError) as exc:
        print(f"validate_mj_tables: {exc}", file=sys.stderr)
        return 1
    print(
        f"MJ_RANDOM_RECORD_VALIDATION PASS records={len(samples)} "
        f"files={len(FAMILIES)} invariant_lengths={len(lengths)} "
        f"class_numbers={len(class_numbers)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
