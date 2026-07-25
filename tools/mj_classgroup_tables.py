#!/usr/bin/env python3
"""Stream Mosunov--Jacobson imaginary-quadratic class-group tables."""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


BLOCK_SIZE = 1 << 28
FILENAME_RE = re.compile(
    r"^(cl(?P<residue>3|7)mod(?P<modulus>8)|"
    r"cl(?P<even_residue>4|8)mod(?P<even_modulus>16))"
    r"\.(?P<block>[0-9]+)\.(?:gz|gzip)$"
)
TSV_HEADER = (
    "D_K\tclass_group_invariants\tclass_number\tp\tp_rank\t"
    "source_filename\tsource_block\tsource_record"
)


class TableFormatError(ValueError):
    """A raw table record or filename does not follow the documented schema."""


@dataclass(frozen=True)
class Source:
    path: Path
    residue: int
    modulus: int
    block: int


@dataclass(frozen=True)
class Record:
    discriminant: int
    class_number: int
    invariants: tuple[int, ...]
    source: Source
    record_number: int

    def p_rank(self, prime: int) -> int:
        return sum(c % prime == 0 for c in self.invariants)


def source_from_path(path: str | Path) -> Source:
    table_path = Path(path)
    match = FILENAME_RE.fullmatch(table_path.name)
    if match is None:
        raise TableFormatError(
            f"{table_path}: expected cl3mod8, cl7mod8, cl4mod16, or "
            "cl8mod16 followed by .BLOCK.gz"
        )
    residue_text = match.group("residue") or match.group("even_residue")
    modulus_text = match.group("modulus") or match.group("even_modulus")
    block = int(match.group("block"))
    if block > 4095:
        raise TableFormatError(f"{table_path}: block must be in 0..4095")
    return Source(table_path, int(residue_text), int(modulus_text), block)


def decode_stream(stream: TextIO, source: Source) -> Iterator[Record]:
    """Decode one decompressed table stream using exact integer arithmetic."""
    discriminant = -(source.block * BLOCK_SIZE + source.residue)
    lower = source.block * BLOCK_SIZE
    upper = (source.block + 1) * BLOCK_SIZE

    for record_number, raw_line in enumerate(stream, 1):
        line = raw_line.strip()
        if not line:
            raise TableFormatError(
                f"{source.path}:{record_number}: blank record"
            )
        fields = line.split()
        if len(fields) < 3:
            raise TableFormatError(
                f"{source.path}:{record_number}: expected a, h, and invariants"
            )
        try:
            values = [int(field, 10) for field in fields]
        except ValueError as exc:
            raise TableFormatError(
                f"{source.path}:{record_number}: non-integer field"
            ) from exc

        delta, class_number, *invariants = values
        if delta < 0:
            raise TableFormatError(
                f"{source.path}:{record_number}: negative discriminant delta"
            )
        if record_number > 1 and delta == 0:
            raise TableFormatError(
                f"{source.path}:{record_number}: repeated discriminant"
            )
        if class_number <= 0 or any(c <= 0 for c in invariants):
            raise TableFormatError(
                f"{source.path}:{record_number}: nonpositive class data"
            )
        if math.prod(invariants) != class_number:
            raise TableFormatError(
                f"{source.path}:{record_number}: invariant product does not "
                f"equal class number {class_number}"
            )

        discriminant -= source.modulus * delta
        magnitude = -discriminant
        if magnitude < lower or magnitude >= upper:
            raise TableFormatError(
                f"{source.path}:{record_number}: |D|={magnitude} outside "
                f"block [{lower},{upper})"
            )
        if magnitude % source.modulus != source.residue:
            raise TableFormatError(
                f"{source.path}:{record_number}: discriminant residue mismatch"
            )
        yield Record(
            discriminant,
            class_number,
            tuple(invariants),
            source,
            record_number,
        )


def read_table(path: str | Path) -> Iterator[Record]:
    source = source_from_path(path)
    try:
        with gzip.open(source.path, mode="rt", encoding="ascii", newline="") as stream:
            yield from decode_stream(stream, source)
    except (gzip.BadGzipFile, EOFError) as exc:
        raise TableFormatError(f"{source.path}: invalid or truncated gzip data") from exc


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1 if divisor == 2 else 2
    return True


def format_record(record: Record, prime: int) -> str:
    invariants = "[" + ",".join(str(c) for c in record.invariants) + "]"
    return "\t".join(
        (
            str(record.discriminant),
            invariants,
            str(record.class_number),
            str(prime),
            str(record.p_rank(prime)),
            record.source.path.name,
            str(record.source.block),
            str(record.record_number),
        )
    )


def filtered_records(
    paths: Iterable[str | Path],
    prime: int,
    min_rank: int,
    min_abs_disc: int,
    max_abs_disc: int,
) -> list[Record]:
    matches: list[Record] = []
    for path in paths:
        for record in read_table(path):
            magnitude = -record.discriminant
            if magnitude > max_abs_disc:
                break
            if magnitude >= min_abs_disc and record.p_rank(prime) >= min_rank:
                matches.append(record)
    matches.sort(key=lambda record: (-record.discriminant, record.source.path.name))
    return matches


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream and filter Mosunov--Jacobson class-group tables."
    )
    parser.add_argument("tables", nargs="+", help="raw .gz/.gzip table files")
    parser.add_argument("--prime", type=int, default=5)
    parser.add_argument("--min-p-rank", type=int, default=3)
    parser.add_argument("--min-abs-disc", type=int, default=1)
    parser.add_argument("--max-abs-disc", type=int, default=(1 << 40) - 1)
    parser.add_argument("--output", help="TSV output path (default: stdout)")
    args = parser.parse_args(argv)
    if not is_prime(args.prime):
        parser.error("--prime must be prime")
    if args.min_p_rank < 0:
        parser.error("--min-p-rank must be nonnegative")
    if args.min_abs_disc < 1:
        parser.error("--min-abs-disc must be positive")
    if args.max_abs_disc < args.min_abs_disc:
        parser.error("--max-abs-disc must not precede --min-abs-disc")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        matches = filtered_records(
            args.tables,
            args.prime,
            args.min_p_rank,
            args.min_abs_disc,
            args.max_abs_disc,
        )
    except (OSError, UnicodeError, TableFormatError) as exc:
        print(f"mj_classgroup_tables: {exc}", file=sys.stderr)
        return 1

    output: TextIO
    try:
        if args.output:
            output = open(args.output, "w", encoding="ascii", newline="\n")
        else:
            output = sys.stdout
    except OSError as exc:
        print(f"mj_classgroup_tables: {exc}", file=sys.stderr)
        return 1
    try:
        print(TSV_HEADER, file=output)
        for record in matches:
            print(format_record(record, args.prime), file=output)
    finally:
        if output is not sys.stdout:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
