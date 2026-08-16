# Transverse rank-one certificates

Finite-algebra certificates for the transverse rank-one mildness
criterion over the full census: for each of the 204 fields of 5-class
rank three with \(|D_K| < 2^{30}\), an explicit finite extension
\(k/\mathbf F_5\) and a point \(x\in\mathbf P^2(k)\) with

```text
rank(D_x) = 1   and   det(B_x) != 0,
```

where \(B_x\) is the 2x2 normal-map matrix built from the polarized
secondary-norm family.  151 fields carry a rational certificate
(k = F_5), 52 require an extension of degree between 2 and 6, and the
single remaining field, D = -781922404, admits no transverse point over
any extension: its norm-degeneracy scheme is one non-reduced closed
point of degree two, recorded here with det B_x = 0 and Jacobian
tangent dimension 1.

`certificates.json` records, per field, every rational rank-one point
and every non-rational closed point of the norm-degeneracy scheme, with kernel
and image bases, the chosen quotient representatives and annihilating
functional, the exact matrix \(B_x\), \(\det B_x\), and the Jacobian
tangent-space dimension used as an independent cross-check
(`det B_x != 0` is equivalent to tangent dimension 0 at every certified
point).  `report.txt` is the human-readable form.

Regenerate deterministically from the committed arithmetic data
(requires Singular for the chart-wise radical/eliminant computation and
the gp of the patched PARI build for reading the certificates):

```sh
python3 tools/transverse_rank_one.py
```

The only inputs are the six verified secondary-norm matrices of each
field, read from its arithmetic certificate in `certificates/p5/`; for the
61 fields with \(|D_K| < 2^{28}\) the committed `result.gp` records
(and, for the principal example,
`certificates/p5/K-2800905-p5/secondary-norms.gp`) are read as well and
must agree with the certificate matrix by matrix.  The expensive
number-field search is never repeated.

`sigma-dimension.json` records the Krull dimension of the minor ideal
of the matrix family over `F_5` for every field, computed by
`tools/sigma_d_dimension.py --p5` from the same certificate export.
All 204 values are 1: every norm-degeneracy scheme is
zero-dimensional, with finitely many closed points.  The chart-wise
radical computation of `certificates.json` is the stronger record
(it lists the closed points themselves); the dimension record is the
uniform statement shared with the `p = 3` and `p = 7` collections.
`sigma-dimension-crosscheck-singular.json` recomputes every
dimension with Singular (`dim(groebner(I))` over `F_5`) and
documents the agreement of the two engines on all 204 fields; the
two share the construction of the minor generators, so the
crosscheck covers the Groebner and dimension step.
