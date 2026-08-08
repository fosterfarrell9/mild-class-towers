#!/usr/bin/env python3
"""Strong-freeness verdicts for the whole block from the verified tensors.

For every field of the block, the truncated two-sided Groebner basis of
the cubic relation ideal is completed to a degree bound and the normal
word counts are compared with the strongly free series 1/(1-3z+3z^3).
Three outcomes are possible:

  STRONGLY_FREE     the completion terminates by the diamond-lemma bound
                    (every basis element of degree at most m and the
                    bound at least 2m-1) and the series agrees through
                    2(s+3)+3, with s the number of automaton states;
  NOT_STRONGLY_FREE a coefficient deviates.  By the coefficientwise
                    inequality for the two series, with equality exactly
                    for a strongly free sequence, this is a proof and
                    not a failed search;
  UNDECIDED         no deviation through the bound, no termination.

The input is the tensor of the C verification, so no arithmetic is
recomputed here: this driver is finite linear algebra over F_3 on data
that ``verifier/verify_all.py`` has already certified.

Deterministic: fixed order xyz, no randomness, no seeds.  Results go to
``results/strong-freeness-block-<n>.json``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)]

import strong_freeness as gb  # noqa: E402

RESULTS = HERE.parent / "results"
TOOLS = HERE.parents[2] / "tools"
BLOCK_SIZE = 2497
SINGULAR = os.environ.get("SINGULAR", "Singular")


def singular_module():
    """Load the Letterplace driver and switch it to characteristic three.

    The driver keeps its own instance of strong_freeness_gb, so the
    characteristic must be set on that instance; setting it on a
    separately imported copy would leave Singular computing over F_5
    while the coefficients still look plausible.
    """
    spec = importlib.util.spec_from_file_location(
        "strong_freeness_singular", TOOLS / "strong_freeness_singular.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.gb.P = 3
    assert module.gb.P == 3
    return module


def latest_verification() -> Path:
    paths = sorted(RESULTS.glob("verification-*.json"))
    for path in reversed(paths):
        data = json.loads(path.read_text())
        if (data.get("all_verified")
                and data.get("certificates_verified") == BLOCK_SIZE):
            return path
    raise FileNotFoundError(
        "no complete block verification in results/; run verifier/verify_all.py")


def next_result() -> Path:
    number = 1
    while (RESULTS / f"strong-freeness-block-{number:03d}.json").exists():
        number += 1
    return RESULTS / f"strong-freeness-block-{number:03d}.json"


class Deviation(Exception):
    def __init__(self, degree: int):
        super().__init__(f"deviation in degree {degree}")
        self.degree = degree


def analyse(item):
    """Verdict for one field.  Runs in a worker process."""
    discriminant, tensor, bound = item
    gb.set_prime(3)
    key = gb.make_order_key("xyz")
    started = time.monotonic()

    def progress(degree, basis):
        series, _ = gb.normal_word_counts(list(basis), degree)
        target = gb.target_series(degree)
        first = next((i for i, pair in enumerate(zip(series, target))
                      if pair[0] != pair[1]), None)
        if first is not None:
            raise Deviation(first)

    deviation = None
    basis = None
    processed = 3
    terminated = False
    try:
        basis, processed, terminated = gb.complete(
            gb.rows_to_polynomials(tensor), key, bound, progress=progress)
    except Deviation as exc:
        deviation = exc.degree

    record = {
        "discriminant": discriminant,
        "degree_bound": bound,
        "first_deviation_degree": deviation,
        "processed_degree": processed,
        "terminated_completion": terminated,
    }
    if deviation is not None:
        record["verdict"] = "NOT_STRONGLY_FREE"
    else:
        leads = list(basis)
        series, states = gb.normal_word_counts(leads, processed)
        target = gb.target_series(processed)
        first = next((i for i, pair in enumerate(zip(series, target))
                      if pair[0] != pair[1]), None)
        record["automaton_states"] = states
        record["basis_elements"] = len(leads)
        record["largest_leading_word_degree"] = max(map(len, leads))
        if first is not None:
            record["first_deviation_degree"] = first
            record["verdict"] = "NOT_STRONGLY_FREE"
        elif terminated:
            probe = 2 * (states + 3) + 3
            long_series, _ = gb.normal_word_counts(leads, probe)
            if long_series == gb.target_series(probe):
                record["verdict"] = "STRONGLY_FREE"
                record["positive_certificate_bound"] = probe
                record["leading_words"] = sorted(leads)
            else:
                record["verdict"] = "NOT_STRONGLY_FREE"
                record["first_deviation_degree"] = next(
                    i for i, pair in enumerate(
                        zip(long_series, gb.target_series(probe)))
                    if pair[0] != pair[1])
        else:
            record["verdict"] = "UNDECIDED"
    record["seconds"] = time.monotonic() - started
    return record


def analyse_singular(item):
    """Same verdict, computed by Singular's Letterplace subsystem."""
    discriminant, tensor, bound = item
    sing = singular_module()
    started = time.monotonic()
    script = sing.build_script(tensor, "xyz", bound)
    completed = subprocess.run([SINGULAR, "-q"], input=script, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    leads = sing.parse_leads(completed.stdout)
    largest = max(len(word) for word in leads)
    terminated = bound >= 2 * largest - 1

    record = {
        "discriminant": discriminant,
        "degree_bound": bound,
        "terminated_completion": terminated,
        "basis_elements": len(leads),
        "largest_leading_word_degree": largest,
    }
    if terminated:
        _, states = gb.normal_word_counts(leads, 1)
        probe = 2 * (states + 3) + 3
        series, states = gb.normal_word_counts(leads, probe)
        target = gb.target_series(probe)
        record["positive_certificate_bound"] = probe
    else:
        series, states = gb.normal_word_counts(leads, bound)
        target = gb.target_series(bound)
    record["automaton_states"] = states
    first = next((n for n in range(len(series)) if series[n] != target[n]),
                 None)
    record["first_deviation_degree"] = first
    if first is not None:
        record["verdict"] = "NOT_STRONGLY_FREE"
    elif terminated:
        record["verdict"] = "STRONGLY_FREE"
        record["leading_words"] = sorted(leads)
    else:
        record["verdict"] = "UNDECIDED"
    record["seconds"] = time.monotonic() - started
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree-bound", type=int, default=13)
    parser.add_argument("--verification", type=Path)
    parser.add_argument("--only", type=int, nargs="*",
                        help="restrict to these discriminants")
    parser.add_argument("--recheck", type=Path,
                        help="restrict to the decided fields of an earlier "
                             "result file and compare the verdicts")
    parser.add_argument("--engine", choices=("python", "singular"),
                        default="python",
                        help="python completes the basis here; singular "
                             "calls Letterplace and reuses the audited "
                             "counting.  Every published verdict should be "
                             "produced by both.")
    parser.add_argument("--workers", type=int,
                        default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    path = (args.verification or latest_verification()).resolve()
    verification = json.loads(path.read_text())
    fields = [(r["discriminant"], r["tensor_3_by_27"])
              for r in verification["records"]]
    claimed = {}
    if args.recheck:
        earlier = json.loads(args.recheck.read_text())
        claimed = {r["discriminant"]: r["verdict"] for r in earlier["records"]
                   if r["verdict"] != "UNDECIDED"}
    wanted = set(args.only or []) | set(claimed)
    if wanted:
        fields = [f for f in fields if f[0] in wanted]
    if args.engine == "singular" and shutil.which(SINGULAR) is None:
        raise SystemExit(f"Singular executable not found: {SINGULAR}")
    print(f"verification={path.name} fields={len(fields)} "
          f"bound={args.degree_bound} engine={args.engine} "
          f"workers={args.workers}", flush=True)

    started = time.monotonic()
    records = []
    counts = {"STRONGLY_FREE": 0, "NOT_STRONGLY_FREE": 0, "UNDECIDED": 0}
    work = [(d, t, args.degree_bound) for d, t in fields]
    worker = analyse_singular if args.engine == "singular" else analyse
    disagreements = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for done, record in enumerate(pool.map(worker, work, chunksize=4), 1):
            records.append(record)
            counts[record["verdict"]] += 1
            expected = claimed.get(record["discriminant"])
            if expected is not None and expected != record["verdict"]:
                disagreements += 1
                print(f"D={record['discriminant']} DISAGREES: "
                      f"{record['verdict']} against {expected}", flush=True)
            if record["verdict"] != "UNDECIDED":
                print(f"D={record['discriminant']:<12} {record['verdict']}"
                      f" deviation={record['first_deviation_degree']}"
                      f" terminated={record['terminated_completion']}",
                      flush=True)
            if done % 200 == 0:
                print(f"... {done}/{len(fields)} {counts}", flush=True)
    elapsed = time.monotonic() - started

    records.sort(key=lambda r: r["discriminant"], reverse=True)
    output = args.output or next_result()
    output.write_text(json.dumps({
        "prime": 3,
        "order": "xyz",
        "degree_bound": args.degree_bound,
        "engine": args.engine,
        "verification": path.name,
        "fields": len(records),
        "randomness_used": False,
        "seeds": [],
        "counts": counts,
        "seconds": elapsed,
        "records": records,
    }, indent=2) + "\n")
    for verdict, number in counts.items():
        print(f"{verdict}: {number}")
    if claimed:
        print(f"agreement with {args.recheck.name}: "
              f"{len(claimed) - disagreements}/{len(claimed)}")
    print(f"seconds={elapsed:.0f}")
    try:
        shown = output.relative_to(HERE.parents[2])
    except ValueError:
        shown = output
    print(f"RESULT={shown}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
