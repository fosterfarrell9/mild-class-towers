#!/usr/bin/env python3
"""Cheap orchestration tests for the resumable audited batch runner."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "run_mildness_batch.py"
SPEC = importlib.util.spec_from_file_location("run_mildness_batch", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
batch = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = batch
SPEC.loader.exec_module(batch)


state = batch.initial_state()
assert len(state["fields"]) == 10
assert [item["D"] for item in state["fields"]] == [
    field.discriminant for field in batch.FIELDS
]
assert state["fields"][0]["polynomial"] == "s^2-s+6990160"
assert state["fields"][4]["polynomial"] == "s^2+13579778"
assert state["fields"][7]["p_primary"] == [125, 5, 5]
assert all(item["state"] == "PENDING" for item in state["fields"])

# Atomic JSON writes leave a complete parseable target and no temporary file.
with tempfile.TemporaryDirectory() as directory_name:
    directory = Path(directory_name)
    target = directory / "state.json"
    batch.atomic_write_json(target, {"generation": 1})
    batch.atomic_write_json(target, {"generation": 2})
    assert json.loads(target.read_text(encoding="utf-8")) == {"generation": 2}
    assert list(directory.glob(".state.json.tmp-*")) == []

# A merged stdout/stderr stream is forwarded live, fully logged, and receives
# heartbeats even while the harmless child is silent.
with tempfile.TemporaryDirectory() as directory_name:
    run_log = Path(directory_name) / "run.log"
    seen: list[str] = []
    heartbeat_count = 0

    def heartbeat() -> str:
        global heartbeat_count
        heartbeat_count += 1
        return f"TEST HEARTBEAT {heartbeat_count}"

    command = [
        sys.executable,
        "-u",
        "-c",
        (
            "import sys,time;"
            "print('dummy stdout',flush=True);"
            "print('dummy stderr',file=sys.stderr,flush=True);"
            "time.sleep(.12);"
            "print('dummy done',flush=True)"
        ),
    ]
    terminal = io.StringIO()
    with contextlib.redirect_stdout(terminal):
        status = batch.run_streamed_process(
            command,
            run_log,
            0.03,
            heartbeat,
            seen.append,
        )
    assert status == 0
    transcript = run_log.read_text(encoding="utf-8")
    assert transcript == "dummy stdout\ndummy stderr\ndummy done\n"
    assert seen == ["dummy stdout", "dummy stderr", "dummy done"]
    assert "dummy stdout" in terminal.getvalue()
    assert "dummy stderr" in terminal.getvalue()
    assert heartbeat_count >= 2
    assert "TEST HEARTBEAT" in terminal.getvalue()

# Interrupted work becomes pending, while completed and failed entries remain
# terminal and therefore are skipped by the ordinary resume work list.
state = batch.initial_state()
state["fields"][0]["state"] = "COMPLETED"
state["fields"][1]["state"] = "RUNNING"
state["fields"][2]["state"] = "FAILED"
batch.prepare_resume(state)
assert state["fields"][0]["state"] == "COMPLETED"
assert state["fields"][1]["state"] == "PENDING"
assert "interrupted_at" in state["fields"][1]
assert state["fields"][2]["state"] == "FAILED"
assert [
    item["D"] for item in state["fields"] if item["state"] == "PENDING"
][0] == batch.FIELDS[1].discriminant

# A field-level failure is persisted with knowledge-valued mildness semantics,
# and the following field remains available for sequential continuation.
with tempfile.TemporaryDirectory() as directory_name:
    directory = Path(directory_name)
    state = batch.initial_state()
    item = state["fields"][0]
    item["state"] = "RUNNING"
    logger = batch.BatchLogger(directory / batch.LOG_NAME)
    with contextlib.redirect_stdout(io.StringIO()):
        batch.mark_failed(
            directory,
            state,
            item,
            time.monotonic(),
            RuntimeError("dummy arithmetic failure"),
            logger,
        )
    assert item["state"] == "FAILED"
    assert item["result_status"] == "ARITHMETIC_COMPUTATION_FAILED"
    assert item["MILD"] == "UNKNOWN"
    assert item["CD"] == "UNKNOWN"
    assert state["fields"][1]["state"] == "PENDING"
    persisted = json.loads(
        (directory / batch.STATE_NAME).read_text(encoding="utf-8")
    )
    assert persisted["fields"][0]["state"] == "FAILED"
    assert persisted["fields"][1]["state"] == "PENDING"

print("MILDNESS_BATCH_RUNNER_TEST PASS")
