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
24  -145367147   (batch; computed by the per-character parallel
                  drivers in parallelization/)
45  -207666763   (batch; likewise)
```

The first nine are exactly the nine smallest absolute discriminants;
beyond rank nine the batch list was not ordered by absolute
discriminant.  The one remaining batch field, D = -109909943 (rank 17,
class number 10000), currently exceeds the memory budget of the
relative-field arithmetic (the PARI stack overflows even at 8 GiB per
character process).  Every computed field has an arithmetic certificate
under `certificate/`.
