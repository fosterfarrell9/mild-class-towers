# Computed examples for p = 5

Layout of this directory:

- `K-18397407/` — the second-smallest example, computed individually
  before the batch infrastructure existed; contains its `result.gp` and
  `run.log`.
- `batch-block0-01/` — the production batch of ten further fields.  One
  subdirectory `D<disc>` per computed field with `result.gp` and
  `run.log`; `batch-state.json` (resumable state, including SHA-256
  hashes of every committed result), `batch.log` (audit trail), and
  `batch.tsv` (summary table) are the provenance record of the runs.
  The batch position numbering used in the log (`FIELD N/10`) is *not*
  the ordering by absolute discriminant.
- `further/` — fields beyond the original batch, computed with the
  per-character parallel drivers in `parallelization/`; one
  subdirectory `D<disc>` with `result.gp` and `run.log` each.
- `transverse-rank-one/` — transverse rank-one certificates for all
  computed fields (see its README).

The principal example \(K=\mathbf Q(\sqrt{-2800905})\),
\(D_K=-11203620\), predates the result-record format; its arithmetic is
documented and certified in `certificate/K-2800905-p5/` (with the six
secondary-norm matrices exported to `secondary-norms.gp` there), and it
is treated in detail in the paper.

Ordered by absolute discriminant (rank within the 61 fields of 5-class
rank three below 2^28), the computed examples are

```text
 1  -11203620    (principal; certificate/K-2800905-p5/)
 2  -18397407    (K-18397407/)
 3  -27960639    (batch)
 4  -35663739    (batch)
 5  -51213139    (batch)
 6  -54319112    (batch)
 7  -61040707    (batch)
 8  -65818135    (batch)
 9  -75949255    (batch)
11  -89017304    (further/)
12  -89218664    (further/)
13  -90903207    (further/)
14  -93121640    (further/)
15  -104545864   (further/)
18  -123482119   (further/)
19  -123779560   (further/)
20  -126740891   (further/; no F_5 Anick witness -- mild by the
                  transverse rank-one criterion, degree-6 point)
24  -145367147   (batch; computed by the per-character parallel
                  drivers in parallelization/)
45  -207666763   (batch; likewise)
```

The first nine are exactly the nine smallest absolute discriminants;
beyond rank nine the batch list was not ordered by absolute
discriminant.  Of the twenty smallest, three remain uncomputed for
memory reasons (the PARI stack of a single character process
overflows the listed ceiling): rank 10, D = -77778287 (12 GiB); rank
16, D = -106660295 (untested, expected similar); rank 17,
D = -109909943 (8 GiB, still a pending batch field).  Every computed
field has an arithmetic certificate under `certificate/`.
