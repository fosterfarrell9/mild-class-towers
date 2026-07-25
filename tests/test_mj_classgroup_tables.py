#!/usr/bin/env python3
"""Focused tests for the Mosunov--Jacobson table decoder."""

from __future__ import annotations

import gzip
import importlib.util
import io
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "mj_classgroup_tables.py"
SPEC = importlib.util.spec_from_file_location("mj_classgroup_tables", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
mj = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mj
SPEC.loader.exec_module(mj)


def expect_error(action, fragment: str) -> None:
    try:
        action()
    except mj.TableFormatError as exc:
        assert fragment in str(exc), (fragment, str(exc))
    else:
        raise AssertionError(f"expected TableFormatError containing {fragment!r}")


def records(text: str, name: str):
    source = mj.source_from_path(name)
    return list(mj.decode_stream(io.StringIO(text), source))


# Published Mosunov example: first and subsequent record in cl4mod16.1.
published = records("0 12160 380 4 4 2\n2 4392 2196 2\n", "cl4mod16.1.gz")
assert [r.discriminant for r in published] == [-268435460, -268435492]
assert published[0].invariants == (380, 4, 4, 2)
assert published[1].class_number == 4392

# First records in each relevant residue family.
assert records("0 1 1\n", "cl3mod8.0.gz")[0].discriminant == -3
assert records("0 1 1\n", "cl7mod8.0.gz")[0].discriminant == -7
assert records("0 1 1\n", "cl4mod16.0.gz")[0].discriminant == -4
assert records("0 1 1\n", "cl8mod16.0.gz")[0].discriminant == -8

# A block transition resets the accumulator at the next file's documented base.
last = records(f"{mj.BLOCK_SIZE // 16 - 1} 1 1\n", "cl4mod16.0.gz")[0]
first = records("0 1 1\n", "cl4mod16.1.gz")[0]
assert last.discriminant == -(mj.BLOCK_SIZE - 12)
assert first.discriminant == -(mj.BLOCK_SIZE + 4)
assert -last.discriminant < mj.BLOCK_SIZE <= -first.discriminant

# The maximum documented block remains safely exact with integer arithmetic.
large = records("0 1 1\n", "cl7mod8.4095.gz")[0]
assert large.discriminant == -(4095 * mj.BLOCK_SIZE + 7)
assert isinstance(large.discriminant, int)

# Rank counts invariant factors divisible by p, including rank > 3.
rank_record = records("0 625 5 5 5 5\n", "cl3mod8.0.gz")[0]
assert rank_record.p_rank(5) == 4

# Malformed, truncated, and inconsistent records fail closed.
expect_error(
    lambda: records("1 5\n", "cl3mod8.0.gz"),
    "expected a, h, and invariants",
)
expect_error(
    lambda: records("x 1 1\n", "cl3mod8.0.gz"),
    "non-integer field",
)
expect_error(
    lambda: records("0 10 5\n", "cl3mod8.0.gz"),
    "invariant product",
)
expect_error(
    lambda: records("-1 1 1\n", "cl3mod8.0.gz"),
    "negative discriminant delta",
)
expect_error(
    lambda: records("0 1 1\n0 1 1\n", "cl3mod8.0.gz"),
    "repeated discriminant",
)
expect_error(
    lambda: records(f"{mj.BLOCK_SIZE // 8} 1 1\n", "cl3mod8.0.gz"),
    "outside block",
)

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "cl3mod8.0.gz"
    with gzip.open(path, "wt", encoding="ascii") as stream:
        stream.write("0 125 5 5 5\n1 1 1\n")
    decoded = list(mj.read_table(path))
    assert decoded[0].p_rank(5) == 3
    assert decoded[0].record_number == 1

    compressed = path.read_bytes()
    path.write_bytes(compressed[:-4])
    expect_error(lambda: list(mj.read_table(path)), "truncated gzip")

print("MJ_CLASSGROUP_TABLE_TEST PASS")
