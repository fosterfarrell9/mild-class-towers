#!/usr/bin/env python3
"""Progressive exact Hilbert comparisons with per-degree timing.

``measured_comparison`` completes the truncated two-sided Groebner basis
of a cubic relation ideal degree by degree, comparing the normal word
counts with the strongly free series 1/(1-3z+3z^3) after every step and
stopping at the first deviation.  It is the shared engine of
``block_witnesses.py``, ``block_sweep.py`` and
``verify_strong_freeness.py``.
"""

from __future__ import annotations

import time
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)]

import strong_freeness as gb  # noqa: E402


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
