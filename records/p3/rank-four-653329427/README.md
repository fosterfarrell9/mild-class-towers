# The rank-four field

One imaginary quadratic field with |D_K| < 2^30 has 3-class rank
exceeding three: D_K = -653329427, with class group of invariants
[210, 3, 3, 3], class number 5670, and elementary 3-class group
(3, 3, 3, 3).  The paper reports it as undecided in Section 5.5;
this directory records the computation behind that report.  The
record sits outside the certificate layer of the paper's
verification appendix: the runs below are exact-arithmetic GP
computations, kept as plain logs.

## Result

The Bockstein matrix of the field has rank 4, so the cone
C_beta = P(ker B) is empty and the p = 3 route admits no candidate
point.  The transversality criterion is therefore not applicable,
no transverse element exists on this route, and the field remains
undecided.

## Contents

- `pol.txt` — the defining polynomial y^2 - y + 163332357.
- `pre.log` — field data: discriminant, class group, 3-rank
  (`pre.gp`).
- `D-matrizen.txt` — the ten secondary norm operators as 4 x 4
  matrices over F_3, labelled by the vector x in the fixed
  character basis: the four basis characters and their six
  pairwise sums, d(d+1)/2 = 10 values, which determine the
  quadratic family x -> D_x (`orakel.gp`).
- `kriterium.log` — the criterion step (`kriterium.gp`): the
  Bockstein matrix B[l, i] = D_i[i, l] has rank 4, hence the cone
  is empty.
- `crosscheck-*.log` — an independent second run of the same
  pipeline from scratch; the ten matrices printed in
  `crosscheck-orakel.log` agree with `D-matrizen.txt` entry for
  entry.
- `pre.gp`, `orakel.gp`, `kriterium.gp` — the scripts that
  produced the logs.  `orakel.gp` is a rank-independent pure-GP
  implementation of the secondary-norm arithmetic of Section 4 of
  the paper; it needs no patched PARI.

The logs are kept verbatim as produced, including their German
labels; the verdict line "KEGEL LEER - Kriterium nicht anwendbar"
reads "cone empty - criterion not applicable".

## Rerunning

```sh
cd records/p3/rank-four-653329427
DISC=-653329427 gp -q pre.gp
POL="y^2 - y + 163332357" P=3 gp -q orakel.gp
gp -q kriterium.gp
```

The oracle writes `D-matrizen.txt` into the working directory, and
`kriterium.gp` reads it from there.  The full run takes a few
minutes on ordinary hardware.
