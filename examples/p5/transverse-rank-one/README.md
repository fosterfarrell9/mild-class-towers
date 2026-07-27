# Transverse rank-one certificates

Finite-algebra certificates showing that each of the nine computed
rank-three fields satisfies the transverse rank-one mildness criterion:
an explicit finite extension \(k/\mathbf F_5\) and a point
\(x\in\mathbf P^2(k)\) with

```text
rank(D_x) = 1   and   det(B_x) != 0,
```

where \(B_x\) is the 2x2 normal-map matrix built from the polarized
secondary-norm family.  Seven fields have rational certificates
(k = F_5); the two formerly open fields require extensions:

```text
-35663739   k of degree 6,  F_5[t]/(t^6+t^5+t^4+1)
-51213139   k of degree 4,  F_5[t]/(t^4+t^3+t^2+1)
```

`certificates.json` records, per field, every rational rank-one point
and every non-rational closed point of the rank-drop scheme, with kernel
and image bases, the chosen quotient representatives and annihilating
functional, the exact matrix \(B_x\), \(\det B_x\), and the Jacobian
tangent-space dimension used as an independent cross-check
(`det B_x != 0` is equivalent to tangent dimension 0 at every certified
point).  `report.txt` is the human-readable form.

Regenerate deterministically from the committed arithmetic data
(requires Singular for the chart-wise radical/eliminant computation):

```sh
python3 tools/transverse_rank_one.py
```

The only inputs are the six verified secondary-norm matrices of each
field (from the committed `result.gp` records, and for the principal
example from `certificate/K-2800905-p5/secondary-norms.gp`); the
expensive number-field search is never repeated.
