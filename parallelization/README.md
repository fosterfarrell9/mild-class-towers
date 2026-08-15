# Per-character parallelization

Parallelizes the expensive arithmetic of a single field by computing the
six prescribed characters

```text
a, b, c, a+b, a+c, b+c
```

in six independent processes.  The proven sequential code in `src/` is
**not modified**: the drivers here link against the unmodified object
files of the main build and call the same exported functions
(`my_secondary_norm_operator` with the exact audit enabled, the same
deterministic input preparation `Buchall` / `my_find_p_gens` /
`my_find_units_mod_p` / `my_find_Ja_vect`, base-BNF certification in
every process, and the same finite-algebra finish
`my_triple_massey_word_matrix` / `my_find_strongly_free_witness`).

Files:

- `character_driver.c` — one process = one character; recomputes the
  (cheap, deterministic) base-field data, certifies the base BNF, and
  computes that character's secondary-norm matrix with the exact audit;
  with `MASSEY_CERTIFICATE_EXPORT` set it writes a well-formed partial
  certificate holding this character's three entries.
- `finish_driver.c` — reads the six matrices, rebuilds the quadratic
  family, and performs exactly the finite-field tail of the sequential
  pipeline (quadratic-scaling checks, tensor identities, rank,
  strong-freeness search), writing a result record with the same schema.
- `run_parallel.py` — orchestrator: spawn six drivers, merge the six
  partial certificates (asserting identical headers) into one
  certificate in the canonical entry order, assemble the result, and
  verify the merged certificate with the shared standalone verifier
  against the fresh result (`RESULT_RECORD_MATCH`).  With
  `--verify-against <committed result.gp>` it additionally requires the
  fresh record to reproduce the committed one entry for entry and the
  certificate to verify against the committed record (regression mode).

Build and run, from this directory (main build first, `make
PARI=<prefix>` in the repository root):

```sh
make PARI="$HOME/.local"
python3 run_parallel.py --polynomial 's^2-s+15260177' \
  --workdir work/D-61040707 \
  --verify-against ../records/p5/batch-block0-01/D-61040707/result.gp
```

Wall-clock for a field drops to roughly the time of its slowest single
character (plus the per-process base-field setup, seconds to minutes).
The `work/` outputs are transient; committed artifacts should continue
to go through the established layouts (`records/p5/...`,
`certificates/p5/K-<n>-p5/`).
