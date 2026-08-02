#!/usr/bin/env python3
"""Progressive exact Hilbert comparisons with per-degree timing."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINITE = ROOT / "experiments" / "p3-finite-algebra"
PILOT = ROOT / "experiments" / "p3-arithmetic-pilot"
SWEEP = ROOT / "experiments" / "p3-delta-sweep"
sys.path[:0] = [str(FINITE), str(PILOT), str(SWEEP)]

import strong_freeness as gb  # noqa: E402
from analyze import MINIMAL_POINTS  # noqa: E402
from analyze_sweep import load_values  # noqa: E402
from reconstruction import reconstruct_from_points  # noqa: E402

RESULTS = HERE / "results"


class ExactDeviation(Exception):
    def __init__(self, degree, series, target):
        super().__init__(f"exact Hilbert deviation in degree {degree}")
        self.degree = degree
        self.series = series
        self.target = target


def measured_comparison(tensor, max_degree: int, label: str):
    gb.set_prime(3)
    key = gb.make_order_key("xyz")
    started = time.monotonic()
    previous = started
    records = []

    def progress(degree, basis):
        nonlocal previous
        now = time.monotonic()
        series, states = gb.normal_word_counts(list(basis), degree)
        target = gb.target_series(degree)
        deviation = next((i for i, values in enumerate(zip(series, target))
                          if values[0] != values[1]), None)
        record = {
            "degree": degree,
            "hilbert_dimension": series[degree],
            "target_dimension": target[degree],
            "matches_through_degree": deviation is None,
            "first_deviation_degree": deviation,
            "basis_elements": len(basis),
            "largest_leading_word_degree": max(map(len, basis)),
            "automaton_states": states,
            "incremental_seconds": now - previous,
            "cumulative_seconds": now - started,
        }
        records.append(record)
        previous = now
        print(f"HILBERT {label} degree={degree} dim={series[degree]} "
              f"target={target[degree]} basis={len(basis)} "
              f"step_s={record['incremental_seconds']:.3f} "
              f"total_s={record['cumulative_seconds']:.3f}", flush=True)
        if deviation is not None:
            raise ExactDeviation(deviation, series, target)

    deviation = None
    terminated = False
    processed = 3
    basis = None
    try:
        basis, processed, terminated = gb.complete(
            gb.rows_to_polynomials(tensor), key, max_degree, progress=progress)
    except ExactDeviation as exc:
        deviation = exc.degree

    elapsed = time.monotonic() - started
    if deviation is None:
        if basis is None:
            raise AssertionError("completion returned no basis")
        final_series, states = gb.normal_word_counts(list(basis), processed)
        target = gb.target_series(processed)
        deviation = next((i for i, values in enumerate(zip(final_series, target))
                          if values[0] != values[1]), None)
    else:
        final_series, target, states = None, None, None
    positive_certificate_bound = None
    if deviation is None and terminated:
        positive_certificate_bound = 2 * (states + 3) + 3
        certificate_series, _ = gb.normal_word_counts(
            list(basis), positive_certificate_bound)
        certificate_target = gb.target_series(positive_certificate_bound)
        deviation = next((i for i, values in enumerate(
            zip(certificate_series, certificate_target))
                          if values[0] != values[1]), None)
    verdict = ("NOT_STRONGLY_FREE" if deviation is not None else
               "STRONGLY_FREE" if terminated else "INCONCLUSIVE")
    return {
        "label": label,
        "requested_max_degree": max_degree,
        "processed_degree": records[-1]["degree"] if records else processed,
        "first_exact_deviation_degree": deviation,
        "verdict": verdict,
        "terminated_completion": terminated,
        "positive_certificate_comparison_bound": positive_certificate_bound,
        "per_degree": records,
        "final_series": final_series,
        "target_series": target,
        "final_automaton_states": states,
        "elapsed_seconds": elapsed,
        "basis_independence": (
            "The generated two-sided ideal depends only on the relation span; "
            "generator changes induce graded free-algebra automorphisms. Any "
            "deviation is therefore an exact basis-independent negative."
        ),
    }


def pilot_tensor():
    return json.loads((PILOT / "results" / "tensor.json").read_text())[
        "tensor_3_by_27"]


def sweep_tensors():
    records = json.loads((SWEEP / "results" / "fields.json").read_text())
    output = []
    for record in records:
        discriminant = record["discriminant"]
        field_dir = SWEEP / "results" / str(abs(discriminant))
        _, _, _, by_point, _ = load_values(field_dir / "matrices.tsv")
        tensor = reconstruct_from_points(MINIMAL_POINTS, by_point)
        output.append((discriminant, tensor))
    return output


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("pilot", "sweep"))
    parser.add_argument("--max-degree", type=int, required=True)
    args = parser.parse_args(argv)
    RESULTS.mkdir(exist_ok=True)
    if args.mode == "pilot":
        result = measured_comparison(pilot_tensor(), args.max_degree,
                                     "D=-3640387")
        (RESULTS / "hilbert-pilot.json").write_text(
            json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return

    records = []
    for discriminant, tensor in sweep_tensors():
        result = measured_comparison(tensor, args.max_degree,
                                     f"D={discriminant}")
        records.append({"discriminant": discriminant, **result})
        (RESULTS / f"hilbert-{abs(discriminant)}.json").write_text(
            json.dumps(records[-1], indent=2) + "\n")
    summary = {
        "uniform_max_degree": args.max_degree,
        "fields": records,
        "exact_negatives": [record["discriminant"] for record in records
                            if record["verdict"] == "NOT_STRONGLY_FREE"],
        "all_match_through_uniform_degree": all(
            record["verdict"] != "NOT_STRONGLY_FREE" for record in records),
        "total_elapsed_seconds": sum(record["elapsed_seconds"]
                                     for record in records),
    }
    (RESULTS / "hilbert-sweep.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
