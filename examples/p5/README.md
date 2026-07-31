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
10  -77778287    (further/; computed across two machines, see
                  certificate/K-77778287-p5/README.md)
11  -89017304    (further/)
12  -89218664    (further/)
13  -90903207    (further/)
14  -93121640    (further/)
15  -104545864   (further/)
16  -106660295   (further/)
18  -123482119   (further/)
19  -123779560   (further/)
20  -126740891   (further/; no F_5 Anick witness -- mild by the
                  transverse rank-one criterion, degree-6 point)
24  -145367147   (batch; computed by the per-character parallel
                  drivers in parallelization/)
45  -207666763   (batch; likewise)
```

Ranks 1 to 16 are exactly the sixteen smallest absolute discriminants;
beyond that the list is no longer an initial segment.  One field of the
twenty smallest remains uncomputed: rank 17, D = -109909943, whose
character processes reached 37.7 GiB and were still running after ten
hours at a ceiling of 32 GiB.

Cost is not governed by h_K.  Rank 16 carries the largest class number
of the list (h = 11000) and took eleven minutes at 0.7 GiB per
character; rank 10 with h = 6000 needed 16.75 GiB and nine hours.  Nor
is it governed by the class groups of the L_x: the six characters of
rank 10 peak at 17564236, 17563508, 17563328, 17563264, 17563272 and
17563800 kbytes, a spread of 0.006 %, and the same uniformity holds for
rank 17 -- yet the L_x differ from character to character.  A
computation that does not
depend on the character dominates both the time and the memory, and
`MASSEY_LOG_LEVEL=2` identifies it: the class group of the auxiliary
field K(zeta_5), which PARI builds inside `bnrclassfield` for the Kummer
construction.  Splitting a field across six processes therefore buys
wall time but no headroom.

Every computed field has an arithmetic certificate under `certificate/`.
