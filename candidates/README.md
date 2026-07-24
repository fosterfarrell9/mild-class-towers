# Imaginary quadratic candidate discovery

Candidate discovery is an inexpensive filter, separate from the audited
secondary-norm and mildness computation. A rank-three candidate is not thereby
claimed to be mild.

## Mathematical convention

The scanner visits integers `n` in the inclusive user range and considers
`D_K=-n` only when PARI's exact `unegisfundamental(n)` test says that `D_K` is a
negative fundamental **field discriminant**. Thus `--min-abs-disc` and
`--max-abs-disc` refer to `|D_K|`, not to a quadratic radicand.

For each accepted discriminant, `quadclassunit(D_K)` computes the class number
and cyclic class-group invariants `[n_1,...,n_t]`. The p-class rank is

```text
dim_F_p Cl(K)/p = #{i : p divides n_i}.
```

This quadratic class-group calculation is a discovery filter. A candidate
later passed to `--example-result` is independently certified and exactly
audited by that pipeline.

## Fresh scan

```sh
./build/massey --scan-candidates \
  --prime 5 \
  --rank 3 \
  --min-abs-disc 1 \
  --max-abs-disc 100000 \
  --checkpoint-every 10000 \
  --progress-seconds 10 \
  --output candidates/p5-r3-1-100000.gp
```

The bounds are inclusive. This example is illustrative; no production-sized
scan is run as part of the implementation.

`SCAN_PROGRESS` lines report current `|D_K|`, interval percentage, fundamental
discriminants examined, candidates, elapsed wall time, and average rate.
`--progress-seconds` controls their cadence and defaults to 10 seconds.

## Checkpoint and resume

The scanner processes increasing `|D_K|`. It atomically writes
`OUTPUT.tmp` and renames it over `OUTPUT`, so a partial write is never presented
as the checkpoint. The record distinguishes:

- the last absolute discriminant considered, whether fundamental or not;
- the number of fundamental discriminants actually examined;
- all candidates already found.

Resume with exactly the same mathematical range and filters:

```sh
./build/massey --scan-candidates \
  --prime 5 --rank 3 \
  --min-abs-disc 1 --max-abs-disc 100000 \
  --output candidates/p5-r3-1-100000.gp \
  --resume
```

`--stop-after N` provides a controlled pause after `N` integers considered and
is useful for schedulers and resume tests. The record remains `IN_PROGRESS`.

## GP schema

Format version 1 is a vector of `[name,value]` pairs:

1. `format_version`
2. `status`: `IN_PROGRESS` or `COMPLETE`
3. `p`
4. `requested_p_rank`
5. `min_abs_discriminant`
6. `max_abs_discriminant`
7. `last_abs_discriminant_considered`
8. `fundamental_discriminants_examined`
9. `candidates_found`
10. `candidates`

Each candidate is:

```text
[D_K, p, class_group_invariants, class_number, p_rank, polynomial_string]
```

Candidate entries are deterministically ordered by increasing `|D_K|`.

## Inspect and hand off

GP inspection:

```sh
gp -q
? r=read("candidates/p5-r3-1-100000.gp");
? r[10][2]
```

Print commands suitable for manual handoff:

```sh
./build/massey --candidate-inputs candidates/p5-r3-1-100000.gp
```

For one selected candidate:

```sh
./build/massey --example-result result.gp 5 's^2+2800905'
```

The scanner never invokes this expensive command itself.

## Disjoint ranges and deterministic merging

Disjoint inclusive ranges can be scanned in separate processes. Until a
dedicated merge command is added, merge only completed records with identical
`p` and requested rank and nonoverlapping ranges:

1. concatenate their candidate vectors;
2. sort by increasing `abs(D_K)`;
3. sum `fundamental_discriminants_examined`;
4. set the combined minimum and maximum to the union bounds;
5. verify `candidates_found` equals the merged vector length.

Keep the original interval records as the reproducibility evidence. Overlapping
ranges must be deduplicated by `D_K` rather than blindly concatenated.

## Tests

```sh
sh tests/test_candidate_scanner.sh
```

The test covers exact discriminant recognition, both known rank-three fields,
a non-rank-three rejection, deterministic ordering, GP round-trip parsing, and
byte-identical uninterrupted versus checkpoint/resume results.
