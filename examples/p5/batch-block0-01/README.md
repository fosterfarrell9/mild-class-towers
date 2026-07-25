# Batch block0-01: ten rank-three candidates from the Mosunov--Jacobson block-zero list

This directory holds the resumable production batch over ten imaginary
quadratic fields with 5-class rank three, selected from
`candidates/mosunov-jacobson/p5-rank-ge3-block0.tsv` and driven by
`tools/run_mildness_batch.py`. Per-field artifacts are written to
`D-<abs(D_K)>/`; `batch-state.json`, `batch.tsv`, and `batch.log` record the
batch-level state.

A field is only marked COMPLETED after the full audited pipeline (base BNF
certification, prescribed-character normalization, AC1/AC2 by exact ideal
arithmetic, independent norm-class recovery) and an independent PARI/GP
post-validation of the result record.

## First run (2026-07-25, aborted)

The first run started 2026-07-25T06:38:25Z and aborted at 08:26:57Z with

```text
BATCH ABORTED (GLOBAL SAFETY FAILURE) | PARI/GP validation failed:
user error: final no-witness result is not exhaustive
```

State at abort:

| index | D_K | state | outcome |
|---|---|---|---|
| 1 | -27960639 | COMPLETED | MILD = PROVED, CD = 2, audit PASS (18/18) |
| 2 | -35663739 | interrupted | audit complete, bounded search found no witness, MILD = UNKNOWN |
| 3--10 | | PENDING | not started |

### Cause of the abort

The abort was a defect in the batch runner, not a mathematical failure and
not a defect in the audited pipeline itself. `validate_result` in
`tools/run_mildness_batch.py` enforced the rule "a no-witness result must be
exhaustive", which is correct for *final* per-field results, but the runner
also applied it to the *intermediate* bounded-search result. A bounded
search (candidate limit 250000) that finds no strongly free basis therefore
crashed the batch before the intended escalation to the exhaustive
GL_3(F_5) search could start. Field 1 was unaffected because its bounded
search found a witness.

### Field 1, D_K = -27960639 (complete and valid)

`D-27960639/result.gp` and `D-27960639/run.log` are the untouched output of
the first run. The run log records base-field `bnfcertify`, the prescribed
Artin normalization for all six characters, and all 18 AC1/AC2 checks with
independently recovered norm classes. Result summary: class group
[40,10,10], class number 4000, cubic rank 3, strongly free basis found
within the bounded search, leading words `["ccb","caa","cab"]`,
MILD = PROVED, CD = 2. The result record was re-validated independently
after the abort (schema, discriminant, class group, doubled-character
checks, cubic rank, status consistency: all pass).

SHA-256 (also recorded in `batch-state.json`):

```text
result.gp  2146eccba19ed2499c83440602988f65214a0bf1fc06a4f245b54c4b6020384d
run.log    53fb35af430f25cf662100ee118c068d9de51f1260685957d31d61cee3283beb
```

### Field 2, D_K = -35663739 (partial; preserved copies)

The run for this field completed the entire arithmetic audit (base BNF
certification, all six characters, AC1/AC2 PASS, independent norm-class
MATCH, cubic rank 3) and then ran only the bounded strong-freeness search,
which found no witness among the first 250000 candidates. The correct
status at that point is:

```text
MILD = UNKNOWN, CD = UNKNOWN
```

The exhaustive GL_3(F_5) search (1488000 candidates) was never run because
of the abort. **No conclusion about non-mildness may be drawn from this
file:** failure to find a strongly free basis in a bounded search proves
nothing; the status is UNKNOWN until the exhaustive search has run.

Because a resumed batch deletes the temporary bounded result and rewrites
`run.log` for interrupted fields, verbatim copies were preserved before
resuming:

```text
aborted-20260725-result-bounded.gp  copy of .result-bounded.gp from the first run
aborted-20260725-run.log            copy of run.log from the first run
```

SHA-256:

```text
aborted-20260725-result-bounded.gp  041ee3b61ae21928ca48784611303a8561c85b85c86d5eabe0a7a161a3e48bdd
aborted-20260725-run.log            9f3024bf7da54377e7de7ad97a9c197a8cbc8d4f6fb62e8eba4ef14ce18059f0
```

The secondary-norm matrices and the cubic relation matrix in the preserved
bounded result passed the full arithmetic audit and remain mathematically
usable; a later complete run of this field supersedes this partial record
for every question it answers.

### Stale entries in the batch-level files

`batch-state.json`, `batch.tsv`, and `batch.log` are committed as the
verbatim snapshot at the moment of the abort. In that snapshot field 2 is
still marked `RUNNING`; the abort happened before the state file was
updated. A resumed run rewrites these three files.
