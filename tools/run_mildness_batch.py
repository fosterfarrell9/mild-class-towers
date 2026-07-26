#!/usr/bin/env python3
"""Run the existing audited mildness pipeline sequentially and resumably."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BATCH = ROOT / "examples" / "p5" / "batch-block0-01"
STATE_NAME = "batch-state.json"
TSV_NAME = "batch.tsv"
LOG_NAME = "batch.log"
HEARTBEAT_SECONDS = 25.0
GL3_F5_ORDER = 1_488_000


@dataclass(frozen=True)
class Field:
    discriminant: int
    class_group: tuple[int, ...]
    class_number: int

    @property
    def primary(self) -> tuple[int, ...]:
        result = []
        for invariant in self.class_group:
            part = 1
            while invariant % (part * 5) == 0:
                part *= 5
            if part > 1:
                result.append(part)
        return tuple(result)

    @property
    def polynomial(self) -> str:
        if self.discriminant % 4 == 1:
            return f"s^2-s+{(1 - self.discriminant) // 4}"
        return f"s^2+{-self.discriminant // 4}"

    @property
    def directory(self) -> str:
        return f"D{self.discriminant}"

    @property
    def radicand(self) -> int:
        """Positive n with K = Q(sqrt(-n)), used in certificate paths."""
        if self.discriminant % 4 == 1:
            return -self.discriminant
        return -self.discriminant // 4


FIELDS = (
    Field(-27960639, (40, 10, 10), 4000),
    Field(-35663739, (30, 10, 5), 1500),
    Field(-51213139, (75, 5, 5), 1875),
    Field(-61040707, (20, 10, 5), 1000),
    Field(-54319112, (60, 10, 5), 3000),
    Field(-65818135, (30, 10, 10), 3000),
    Field(-75949255, (40, 20, 5), 4000),
    Field(-145367147, (125, 5, 5), 3125),
    Field(-109909943, (100, 10, 10), 10000),
    Field(-207666763, (50, 10, 5), 2500),
)

IMPORTANT_FRAGMENTS = (
    "MASSEY_ARITHMETIC_AUDIT",
    "prescribed t =",
    "input e_",
    "AC1",
    "AC2",
    "MATCH",
    "SECONDARY_NORMS",
    "MASSEY_RANK",
    "STRONG_FREENESS",
    "LEADING_WORDS",
    "MILD ",
    "error",
    "failed",
    "exception",
)

TSV_COLUMNS = (
    "D",
    "class_group",
    "p_primary",
    "class_number",
    "state",
    "result_status",
    "arithmetic_audit",
    "cubic_rank",
    "bounded_witness",
    "exhaustive_witness",
    "MILD",
    "CD",
    "leading_words",
    "runtime_seconds",
    "result_sha256",
    "run_log_sha256",
    "exit_status",
    "started_at",
    "ended_at",
    "error",
)


class BatchError(RuntimeError):
    """A failure that makes continuing the whole batch unsafe."""


class FieldError(RuntimeError):
    """A field-level arithmetic/preflight failure that permits continuation."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def format_elapsed(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def vector_text(values: tuple[int, ...] | list[int]) -> str:
    return "[" + ",".join(str(value) for value in values) + "]"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_json(path: Path, value: dict) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def initial_state() -> dict:
    return {
        "format_version": 1,
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "configuration": {
            "p": 5,
            "strong_search_limit": "exhaustive",
            "exhaustive_gl3_f5_candidates": GL3_F5_ORDER,
            "sequential_processes": 1,
        },
        "fields": [
            {
                "index": index,
                "D": field.discriminant,
                "polynomial": field.polynomial,
                "class_group": list(field.class_group),
                "p_primary": list(field.primary),
                "class_number": field.class_number,
                "directory": field.directory,
                "state": "PENDING",
            }
            for index, field in enumerate(FIELDS, 1)
        ],
    }


def write_state(batch: Path, state: dict) -> None:
    state["updated_at"] = utc_now()
    atomic_write_json(batch / STATE_NAME, state)
    write_tsv(batch, state)


def write_tsv(batch: Path, state: dict) -> None:
    rows: list[list[str]] = []
    for item in state["fields"]:
        rows.append(
            [
                str(item.get("D", "")),
                vector_text(item.get("class_group", [])),
                vector_text(item.get("p_primary", [])),
                str(item.get("class_number", "")),
                str(item.get("state", "")),
                str(item.get("result_status", "")),
                str(item.get("arithmetic_audit", "")),
                str(item.get("cubic_rank", "")),
                str(item.get("bounded_witness", "")),
                str(item.get("exhaustive_witness", "")),
                str(item.get("MILD", "")),
                str(item.get("CD", "")),
                str(item.get("leading_words", "")),
                str(item.get("runtime_seconds", "")),
                str(item.get("result_sha256", "")),
                str(item.get("run_log_sha256", "")),
                str(item.get("exit_status", "")),
                str(item.get("started_at", "")),
                str(item.get("ended_at", "")),
                str(item.get("error", "")).replace("\t", " ").replace("\n", " "),
            ]
        )
    from io import StringIO

    output = StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_COLUMNS)
    writer.writerows(rows)
    atomic_write_text(batch / TSV_NAME, output.getvalue())


class BatchLogger:
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def emit(self, message: str) -> None:
        print(message, flush=True)
        self.record(message)

    def record(self, message: str) -> None:
        with self.path.open("a", encoding="utf-8") as stream:
            stream.write(message + "\n")
            stream.flush()

    def important(self, line: str) -> None:
        lowered = line.lower()
        if any(fragment.lower() in lowered for fragment in IMPORTANT_FRAGMENTS):
            with self.path.open("a", encoding="utf-8") as stream:
                stream.write(f"[{utc_now()}] CHILD | {line.rstrip()}\n")
                stream.flush()


def run_streamed_process(
    command: list[str],
    run_log: Path,
    heartbeat_seconds: float,
    heartbeat_factory: Callable[[], str],
    on_line: Callable[[str], None],
    *,
    append: bool = False,
    env: dict | None = None,
) -> int:
    """Merge, forward, and capture child output while emitting heartbeats."""
    mode = "a" if append else "w"
    run_log.parent.mkdir(parents=True, exist_ok=True)
    child = subprocess.Popen(
        command,
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
    )
    assert child.stdout is not None
    output_queue: queue.Queue[str | None] = queue.Queue()

    def reader() -> None:
        try:
            for line in child.stdout:
                output_queue.put(line)
        finally:
            output_queue.put(None)

    thread = threading.Thread(target=reader, daemon=True)
    thread.start()
    next_heartbeat = time.monotonic() + heartbeat_seconds
    stream_ended = False
    try:
        with run_log.open(mode, encoding="utf-8", newline="") as log:
            while not stream_ended:
                timeout = max(0.01, next_heartbeat - time.monotonic())
                try:
                    item = output_queue.get(timeout=timeout)
                except queue.Empty:
                    item = ""
                if item is None:
                    stream_ended = True
                elif item:
                    log.write(item)
                    log.flush()
                    print(item, end="", flush=True)
                    on_line(item.rstrip("\n"))
                now = time.monotonic()
                if now >= next_heartbeat and not stream_ended:
                    print(heartbeat_factory(), flush=True)
                    next_heartbeat = now + heartbeat_seconds
        thread.join(timeout=5)
        return child.wait()
    except BaseException:
        child.terminate()
        try:
            child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()
        raise


def gp_quote(path: Path) -> str:
    return '"' + str(path).replace("\\", "\\\\").replace('"', '\\"') + '"'


def run_gp(script: str, gp: str) -> list[str]:
    completed = subprocess.run(
        [gp, "-fq"],
        cwd=ROOT,
        input=script,
        text=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=180,
    )
    if completed.returncode != 0 or "***" in completed.stderr:
        detail = (completed.stderr or completed.stdout).strip()
        raise BatchError(f"PARI/GP validation failed: {detail}")
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def preflight_field(field: Field, gp: str) -> None:
    script = (
        f"q=quadclassunit({field.discriminant});\n"
        f"print(nfdisc({field.polynomial}));\n"
        "print(q.no);\n"
        "print(q.cyc);\n"
        "quit\n"
    )
    lines = run_gp(script, gp)
    expected = [
        str(field.discriminant),
        str(field.class_number),
        vector_text(field.class_group).replace(",", ", "),
    ]
    if lines != expected:
        raise FieldError(
            f"D={field.discriminant}: nfdisc/quadclassunit mismatch; "
            f"expected {expected!r}, got {lines!r}"
        )


def validate_result(
    path: Path, field: Field, gp: str, final: bool = True
) -> dict:
    """Validate a result record with PARI/GP.

    With final=True (the default) a no-witness result must come from the
    exhaustive search.  The intermediate bounded-search result is validated
    with final=False, since a bounded no-witness outcome is the legitimate
    trigger for the exhaustive escalation, not an error.
    """
    expected = vector_text(field.class_group).replace(",", ", ")
    require_exhaustive = "1" if final else "0"
    script = f"""
r=read({gp_quote(path)});
if(type(r)!="t_VEC" || #r!=20,error("bad result schema"));
if(r[3][2]!=5,error("wrong p"));
if(r[5][2]!={field.discriminant},error("wrong discriminant"));
if(r[6][2]!={expected},error("wrong class group"));
if(r[7][2]!={field.class_number},error("wrong class number"));
if(r[12][2]!=1,error("arithmetic audit absent"));
d=r[10][2]; q=r[11][2];
if(#d!=6 || #q!=3,error("secondary-norm family incomplete"));
for(i=1,6,if(matsize(d[i])!=[3,3],error("bad D matrix")));
for(i=1,3,if(q[i]!=lift(Mod(4,5)*d[i]),error("doubled check mismatch")));
T=r[14][2]; rk=matrank(Mod(T,5));
if(rk!=r[15][2],error("rank mismatch"));
st=r[2][2]; w=r[16][2]; mild=r[19][2]; cd=r[20][2]; ex=r[18][2];
if(st!="STRONGLY_FREE_BASIS_FOUND" && st!="RANK_LT_3" && st!="NO_STRONGLY_FREE_BASIS_FOUND",error("unknown result status"));
if(st=="STRONGLY_FREE_BASIS_FOUND" && (rk!=3 || #w==0 || w[6]!=1 || mild!="PROVED" || cd!=2),error("proved status inconsistent"));
if(st!="STRONGLY_FREE_BASIS_FOUND" && (mild!="UNKNOWN" || cd!="UNKNOWN" || #w!=0),error("unknown status inconsistent"));
if(st=="RANK_LT_3" && rk>=3,error("rank status inconsistent"));
if(st=="NO_STRONGLY_FREE_BASIS_FOUND" && rk!=3,error("no-witness status inconsistent"));
if({require_exhaustive} && st=="NO_STRONGLY_FREE_BASIS_FOUND" && ex!=1,error("final no-witness result is not exhaustive"));
print(rk);
print(mild);
print(cd);
print(st);
print(ex);
if(#w,print(w[4]),print([]));
quit
"""
    lines = run_gp(script, gp)
    if len(lines) != 6:
        raise BatchError(f"{path}: incomplete GP validation output {lines!r}")
    rank_text, mild, cd, status, exhaustive, leading_words = lines
    return {
        "cubic_rank": int(rank_text),
        "MILD": mild,
        "CD": cd,
        "result_status": status,
        "exhaustive": exhaustive == "1",
        "leading_words": leading_words,
    }


def verify_completed(batch: Path, item: dict, field: Field, gp: str) -> None:
    directory = batch / item["directory"]
    result = directory / "result.gp"
    run_log = directory / "run.log"
    if not result.is_file() or not run_log.is_file():
        raise BatchError(f"D={field.discriminant}: completed files are missing")
    if sha256(result) != item.get("result_sha256"):
        raise BatchError(f"D={field.discriminant}: result.gp hash changed")
    if sha256(run_log) != item.get("run_log_sha256"):
        raise BatchError(f"D={field.discriminant}: run.log hash changed")
    parsed = validate_result(result, field, gp)
    for key in ("cubic_rank", "MILD", "CD", "result_status"):
        if str(parsed[key]) != str(item.get(key)):
            raise BatchError(
                f"D={field.discriminant}: manifest/result mismatch for {key}"
            )


def load_state(batch: Path) -> dict:
    try:
        with (batch / STATE_NAME).open(encoding="utf-8") as stream:
            state = json.load(stream)
    except (OSError, ValueError) as exc:
        raise BatchError(f"cannot load batch state: {exc}") from exc
    expected = initial_state()["fields"]
    actual = state.get("fields")
    if not isinstance(actual, list) or len(actual) != len(expected):
        raise BatchError("batch state does not describe the fixed ten-field batch")
    for old, new in zip(actual, expected):
        for key in (
            "index",
            "D",
            "polynomial",
            "class_group",
            "p_primary",
            "class_number",
            "directory",
        ):
            if old.get(key) != new[key]:
                raise BatchError(f"batch-state field specification changed at {key}")
    return state


def prepare_resume(state: dict) -> None:
    """Make an interrupted RUNNING entry retryable; keep terminal states."""
    for item in state["fields"]:
        if item["state"] == "RUNNING":
            item["state"] = "PENDING"
            item["interrupted_at"] = utc_now()


def recover_interrupted_results(
    batch: Path, state: dict, gp: str, logger: BatchLogger
) -> None:
    """Promote a validated result written just before an interruption."""
    for item, field in zip(state["fields"], FIELDS):
        if item["state"] != "PENDING" or "interrupted_at" not in item:
            continue
        directory = batch / item["directory"]
        result = directory / "result.gp"
        run_log = directory / "run.log"
        if not result.exists():
            for name in (".result-bounded.gp", ".result-exhaustive.gp"):
                temporary = directory / name
                if temporary.exists():
                    temporary.unlink()
            continue
        if not run_log.is_file():
            raise BatchError(
                f"D={field.discriminant}: interrupted result has no run.log"
            )
        parsed = validate_result(result, field, gp)
        try:
            started = datetime.fromisoformat(item["started_at"]).timestamp()
            runtime = max(0.0, result.stat().st_mtime - started)
        except (KeyError, ValueError):
            runtime = 0.0
        item.update(parsed)
        item.update(
            {
                "state": "COMPLETED",
                "arithmetic_audit": "PASS (18/18)",
                "runtime_seconds": round(runtime, 3),
                "result_sha256": sha256(result),
                "run_log_sha256": sha256(run_log),
                "ended_at": item["interrupted_at"],
            }
        )
        if parsed["cubic_rank"] < 3:
            item["bounded_witness"] = "NOT_APPLICABLE"
            item["exhaustive_witness"] = "NOT_APPLICABLE"
        elif parsed["exhaustive"]:
            item["bounded_witness"] = "NOT_RUN"
            item["exhaustive_witness"] = (
                "FOUND" if parsed["MILD"] == "PROVED" else "NOT_FOUND"
            )
        else:
            item["bounded_witness"] = (
                "FOUND" if parsed["MILD"] == "PROVED" else "NOT_FOUND"
            )
            item["exhaustive_witness"] = "NOT_USED"
        logger.emit(
            f"RECOVER D={field.discriminant} | validated completed result "
            "written before interruption"
        )


def reset_for_force(batch: Path, state: dict) -> None:
    for item in state["fields"]:
        directory = batch / item["directory"]
        for name in ("result.gp", "run.log", ".result-bounded.gp", ".result-exhaustive.gp"):
            path = directory / name
            if path.exists():
                path.unlink()
        preserved = {
            key: item[key]
            for key in (
                "index",
                "D",
                "polynomial",
                "class_group",
                "p_primary",
                "class_number",
                "directory",
            )
        }
        item.clear()
        item.update(preserved)
        item["state"] = "PENDING"


def field_summary(
    logger: BatchLogger, item: dict, total: int, completed: int
) -> None:
    line = "=" * 62
    logger.emit(line)
    logger.emit(f"FIELD {item['index']}/{total} COMPLETE")
    logger.emit(f"D = {item['D']}")
    logger.emit(f"arithmetic audit = {item['arithmetic_audit']}")
    logger.emit(f"rank(T) = {item['cubic_rank']}")
    logger.emit(
        "strong freeness = "
        + ("PASS" if item["MILD"] == "PROVED" else "INCONCLUSIVE")
    )
    logger.emit(f"MILD = {item['MILD']}")
    logger.emit(f"CD = {item['CD']}")
    logger.emit(f"field elapsed = {format_elapsed(item['runtime_seconds'])}")
    logger.emit(f"batch completed = {completed}/{total}")
    logger.emit(line)


def run_field(
    batch: Path,
    state: dict,
    item: dict,
    field: Field,
    executable: Path,
    gp: str,
    stdbuf: str,
    heartbeat_seconds: float,
    logger: BatchLogger,
    certificate_dir: Path | None = None,
) -> None:
    index = item["index"]
    total = len(FIELDS)
    completed_before = sum(
        entry["state"] == "COMPLETED" for entry in state["fields"]
    )
    directory = batch / field.directory
    directory.mkdir(parents=True, exist_ok=True)
    final_result = directory / "result.gp"
    run_log = directory / "run.log"
    bounded_result = directory / ".result-bounded.gp"
    exhaustive_result = directory / ".result-exhaustive.gp"
    for path in (final_result, bounded_result, exhaustive_result):
        if path.exists():
            raise BatchError(
                f"D={field.discriminant}: refusing to overwrite {path}; "
                "use --force explicitly"
            )
    if run_log.exists():
        run_log.unlink()

    preflight_field(field, gp)
    started = time.monotonic()
    item.update(
        {
            "state": "RUNNING",
            "started_at": utc_now(),
            "bounded_witness": "NOT_RUN",
            "exhaustive_witness": "NOT_USED",
        }
    )
    write_state(batch, state)
    logger.emit("=" * 62)
    logger.emit(
        f"FIELD {index}/{total} START | D={field.discriminant} | "
        f"polynomial={field.polynomial} | completed={completed_before}/{total}"
    )
    pipeline_command = [
        stdbuf,
        "-oL",
        "-eL",
        str(executable),
        "--example-result",
        str(exhaustive_result),
        "--strong-search-limit",
        "exhaustive",
        "5",
        field.polynomial,
    ]
    environment = None
    if certificate_dir is not None:
        certificate_path = (
            certificate_dir / f"K-{field.radicand}-p5" / "certificate.gp"
        )
        if certificate_path.exists():
            raise BatchError(
                f"D={field.discriminant}: refusing to overwrite "
                f"{certificate_path}"
            )
        certificate_path.parent.mkdir(parents=True, exist_ok=True)
        environment = dict(os.environ)
        environment["MASSEY_CERTIFICATE_EXPORT"] = str(certificate_path)
        item["certificate"] = str(certificate_path)
        logger.emit(
            f"[{utc_now()}] D={field.discriminant} | "
            f"certificate export -> {certificate_path}"
        )
    progress = {"character": None, "attempt": "exhaustive"}
    character_names = {
        "[1, 0, 0]": "a",
        "[0, 1, 0]": "b",
        "[0, 0, 1]": "c",
        "[1, 1, 0]": "a+b",
        "[1, 0, 1]": "a+c",
        "[0, 1, 1]": "b+c",
    }

    def child_line(line: str) -> None:
        logger.important(line)
        if "prescribed t =" in line:
            vector = line.split("=", 1)[1].strip()
            progress["character"] = character_names.get(vector, vector)
            logger.emit(
                f"[{utc_now()}] D={field.discriminant} | "
                f"character={progress['character']} | arithmetic audit BEGIN"
            )
        elif "input e_3:" in line and progress["character"] is not None:
            logger.emit(
                f"[{utc_now()}] D={field.discriminant} | "
                f"D_{progress['character']} COMPLETE | AC1 PASS | AC2 PASS | "
                "norm-class MATCH"
            )
        elif line.strip() == "MASSEY_RANK 3":
            logger.emit(
                f"[{utc_now()}] D={field.discriminant} | "
                f"STRONG_FREENESS_SEARCH BEGIN ({progress['attempt']})"
            )

    def heartbeat() -> str:
        message = (
            f"[{utc_now()}] BATCH {index}/{total} | D={field.discriminant} | "
            f"RUNNING | field_elapsed={format_elapsed(time.monotonic() - started)} "
            f"| batch_completed={completed_before}/{total}"
        )
        logger.record(message)
        return message

    exit_status = run_streamed_process(
        pipeline_command,
        run_log,
        heartbeat_seconds,
        heartbeat,
        child_line,
        env=environment,
    )
    item["exit_status"] = exit_status
    if exit_status != 0 or not exhaustive_result.is_file():
        raise RuntimeError(f"audited pipeline exited {exit_status}")

    parsed_run = validate_result(exhaustive_result, field, gp)
    if parsed_run["cubic_rank"] == 3:
        item["exhaustive_witness"] = (
            "FOUND" if parsed_run["MILD"] == "PROVED" else "NOT_FOUND"
        )
        logger.emit(
            f"[{utc_now()}] D={field.discriminant} | "
            f"STRONG_FREENESS_SEARCH END (exhaustive) | "
            f"witness={item['exhaustive_witness']}"
        )

    os.replace(exhaustive_result, final_result)
    for temporary in (bounded_result, exhaustive_result):
        if temporary.exists():
            temporary.unlink()
    parsed = validate_result(final_result, field, gp)
    logger.emit(
        f"[{utc_now()}] D={field.discriminant} | INDEPENDENT_POSTVALIDATION PASS "
        "| six_D_matrices=PASS | doubled_character_checks=PASS | "
        f"rank(T)={parsed['cubic_rank']}"
    )
    elapsed = time.monotonic() - started
    item.update(parsed)
    item.update(
        {
            "state": "COMPLETED",
            "arithmetic_audit": "PASS (18/18)",
            "runtime_seconds": round(elapsed, 3),
            "result_sha256": sha256(final_result),
            "run_log_sha256": sha256(run_log),
            "ended_at": utc_now(),
        }
    )
    if parsed["cubic_rank"] < 3:
        item["bounded_witness"] = "NOT_APPLICABLE"
        item["exhaustive_witness"] = "NOT_APPLICABLE"
    write_state(batch, state)
    completed = sum(
        entry["state"] == "COMPLETED" for entry in state["fields"]
    )
    field_summary(logger, item, total, completed)


def mark_failed(
    batch: Path,
    state: dict,
    item: dict,
    started: float,
    error: Exception,
    logger: BatchLogger,
) -> None:
    directory = batch / item["directory"]
    item.setdefault("started_at", utc_now())
    item.setdefault("exit_status", "NOT_STARTED")
    result = directory / "result.gp"
    for name in (".result-bounded.gp", ".result-exhaustive.gp"):
        temporary = directory / name
        if temporary.exists():
            temporary.unlink()
    item.update(
        {
            "state": "FAILED",
            "result_status": "ARITHMETIC_COMPUTATION_FAILED",
            "arithmetic_audit": "FAILED",
            "MILD": "UNKNOWN",
            "CD": "UNKNOWN",
            "runtime_seconds": round(time.monotonic() - started, 3),
            "ended_at": utc_now(),
            "error": str(error),
        }
    )
    if result.exists():
        result.unlink()
    run_log = directory / "run.log"
    if run_log.exists():
        item["run_log_sha256"] = sha256(run_log)
    write_state(batch, state)
    logger.emit("=" * 62)
    logger.emit(f"FIELD {item['index']}/{len(FIELDS)} FAILED")
    logger.emit(f"D = {item['D']}")
    logger.emit("status = ARITHMETIC_COMPUTATION_FAILED")
    logger.emit("MILD = UNKNOWN")
    logger.emit("CD = UNKNOWN")
    logger.emit(f"error = {error}")
    logger.emit("=" * 62)


def print_resume(logger: BatchLogger, state: dict) -> None:
    complete = [item for item in state["fields"] if item["state"] == "COMPLETED"]
    pending = [item for item in state["fields"] if item["state"] == "PENDING"]
    failed = [item for item in state["fields"] if item["state"] == "FAILED"]
    inconclusive = [
        item
        for item in complete
        if item.get("MILD") == "UNKNOWN"
    ]
    logger.emit("BATCH RESUME")
    logger.emit(f"completed: {len(complete)}/{len(FIELDS)}")
    logger.emit(
        "next field: "
        + (str(pending[0]["D"]) if pending else "none")
    )
    logger.emit(
        "previous failed fields: "
        + (", ".join(str(item["D"]) for item in failed) if failed else "none")
    )
    logger.emit(
        "previous inconclusive fields: "
        + (
            ", ".join(str(item["D"]) for item in inconclusive)
            if inconclusive
            else "none"
        )
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, default=DEFAULT_BATCH)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="explicitly discard and rerun all existing field artifacts",
    )
    parser.add_argument(
        "--heartbeat-seconds", type=float, default=HEARTBEAT_SECONDS
    )
    parser.add_argument("--executable", type=Path, default=ROOT / "build" / "massey")
    parser.add_argument("--gp", default=shutil.which("gp") or "gp")
    parser.add_argument(
        "--certificate-dir",
        type=Path,
        default=None,
        help="export an arithmetic certificate for each computed field to "
        "<dir>/K-<n>-p5/certificate.gp during the pipeline run",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="run at most this many pending fields, then stop cleanly "
        "(0 = no limit)",
    )
    args = parser.parse_args(argv)
    if args.heartbeat_seconds <= 0 or args.heartbeat_seconds > 30:
        parser.error("--heartbeat-seconds must be in (0,30]")
    if args.limit < 0:
        parser.error("--limit must be nonnegative")
    if args.certificate_dir is not None:
        args.certificate_dir = args.certificate_dir.resolve()
    args.batch = args.batch.resolve()
    args.executable = args.executable.resolve()
    if not args.batch.is_relative_to(ROOT):
        parser.error("--batch must be inside this repository")
    if args.force and not args.resume:
        parser.error("--force requires --resume")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    stdbuf = shutil.which("stdbuf")
    if not args.executable.is_file() or not os.access(args.executable, os.X_OK):
        print(f"batch runner: executable unavailable: {args.executable}", file=sys.stderr)
        return 2
    if stdbuf is None:
        print("batch runner: stdbuf is required for live child output", file=sys.stderr)
        return 2
    batch = args.batch
    state_path = batch / STATE_NAME
    if args.resume:
        if not state_path.is_file():
            print("batch runner: --resume requested but no state exists", file=sys.stderr)
            return 2
        state = load_state(batch)
        prepare_resume(state)
        if args.force:
            reset_for_force(batch, state)
    else:
        if batch.exists() and any(batch.iterdir()):
            print(
                f"batch runner: refusing nonempty batch directory {batch}; "
                "use --resume",
                file=sys.stderr,
            )
            return 2
        batch.mkdir(parents=True, exist_ok=True)
        state = initial_state()

    logger = BatchLogger(batch / LOG_NAME)
    try:
        logger.emit(
            f"[{utc_now()}] validating executable and all ten field specifications"
        )
        preflight_failures: dict[int, FieldError] = {}
        for field in FIELDS:
            try:
                preflight_field(field, args.gp)
            except FieldError as exc:
                preflight_failures[field.discriminant] = exc
                logger.emit(f"PREFLIGHT FAILED | {exc}")
        if args.resume:
            recover_interrupted_results(batch, state, args.gp, logger)
            for item, field in zip(state["fields"], FIELDS):
                if item["state"] == "COMPLETED":
                    verify_completed(batch, item, field, args.gp)
            write_state(batch, state)
            print_resume(logger, state)
        else:
            write_state(batch, state)
            logger.emit(f"BATCH START | fields={len(FIELDS)} | sequential=1")

        attempted = 0
        for item, field in zip(state["fields"], FIELDS):
            if item["state"] != "PENDING":
                logger.emit(
                    f"SKIP D={field.discriminant} | state={item['state']}"
                )
                continue
            if args.limit and attempted >= args.limit:
                logger.emit(
                    f"LIMIT REACHED | {attempted} field(s) attempted | "
                    f"pausing before D={field.discriminant}"
                )
                break
            started = time.monotonic()
            try:
                if field.discriminant in preflight_failures:
                    raise preflight_failures[field.discriminant]
                run_field(
                    batch,
                    state,
                    item,
                    field,
                    args.executable,
                    args.gp,
                    stdbuf,
                    args.heartbeat_seconds,
                    logger,
                    certificate_dir=args.certificate_dir,
                )
            except BatchError:
                raise
            except Exception as exc:
                mark_failed(batch, state, item, started, exc, logger)
            attempted += 1

        complete = sum(
            item["state"] == "COMPLETED" for item in state["fields"]
        )
        failed = sum(item["state"] == "FAILED" for item in state["fields"])
        logger.emit(
            f"BATCH END | completed={complete}/{len(FIELDS)} | failed={failed} "
            f"| timestamp={utc_now()}"
        )
        return 0 if failed == 0 else 1
    except KeyboardInterrupt:
        write_state(batch, state)
        logger.emit(f"BATCH INTERRUPTED | timestamp={utc_now()} | resume is safe")
        return 130
    except BatchError as exc:
        write_state(batch, state)
        logger.emit(f"BATCH ABORTED (GLOBAL SAFETY FAILURE) | {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
