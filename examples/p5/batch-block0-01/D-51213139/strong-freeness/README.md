# Strong-freeness evidence for D_K = -51213139

The audited record `../result.gp` establishes cubic rank three, so the
rows of its `cubic_relation_matrix` are the complete cubic initial
relation space R_3 of the tower group.  The exhaustive GL_3(F_5) search
of the audited pipeline found no basis with combinatorially free leading
words, so the mildness status of this field is UNKNOWN.  The files in
this directory record the further, purely finite-algebra evidence quoted
in the paper.

## Files and provenance

`python-deglex-<order>.log` — degree-truncated noncommutative Groebner
(diamond-lemma) completion by the repository tool, one file per
degree-lexicographic letter order:

```sh
python3 tools/strong_freeness_gb.py \
  --result examples/p5/batch-block0-01/D-51213139/result.gp \
  --order <order> --max-degree 12 --json
```

`singular-deglex-<order>.log` — the same computation through the
independent Singular/Letterplace engine:

```sh
python3 tools/strong_freeness_singular.py \
  --result examples/p5/batch-block0-01/D-51213139/result.gp \
  --order <order> --degree-bound 13 --json
```

## What the results show

In every one of the six orders and for both engines the completion does
not terminate within the degree bound, and the Hilbert series of the
cubic relation algebra agrees with the strongly free prediction
1/(1-3z+3z^3) in every verified degree: through degree 12 (Python
engine) and through degree 13 (Singular engine).  The verdict field of
every JSON summary is INCONCLUSIVE_SERIES_MATCHES.

## What may not be concluded

Agreement of the Hilbert series through a finite degree does not prove
strong freeness, and failure of the certificate searches does not prove
non-mildness.  A deviation in any verified degree would have disproved
strong freeness; none occurred.  The mildness status of this field
therefore remains open.
