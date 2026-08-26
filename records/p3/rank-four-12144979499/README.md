# A mild field of 3-class rank four (preliminary)

The imaginary quadratic field with D_K = -12144979499 has class group of
invariants [1278, 3, 3, 3], class number 34506, and elementary quotient of
3-class group (9, 3, 3, 3).  Its Bockstein cone is a single
rational point, and at that point both parts of the transversality
criterion hold:

```
rk D_x = 2 = d-2            (rank condition)
Theta matrix = [1, 1; 0, 0; 0, 1],  rank 2  (surjective)
```

The point is a transverse element over F_3 itself, so by the
transversality criterion of the paper the 3-class tower group of
the field is mild, of cohomological dimension 2 --- subject to the
status below.  This is one of three such fields found on
26 August 2026; see the repository README.

## Status

These are preliminary data.  The transversality verdict is exact
linear algebra over class-group data computed with PARI's bnfinit
under GRH; the run was reproduced from scratch with byte-identical
results, and three independent probes confirm it.  The field has
not yet passed the certification layer that backs every result of
the paper.

## Contents

- `pol.txt`, `pre.log` --- defining polynomial and field data.
- `D-matrizen.txt` --- the ten secondary norm operators as 4 x 4
  matrices over F_3 (the four basis characters and their six
  pairwise sums), computed by `orakel-capped-red.gp`.
- `kriterium.log` --- the criterion run: Bockstein matrix of rank 3,
  the single cone point, the rank condition, and the Theta matrix
  (`kriterium.gp`).
- `crosscheck-*.log` --- an independent second run of the same
  pipeline from scratch; the ten matrices printed in
  `crosscheck-orakel.log` agree with `D-matrizen.txt` entry for
  entry.
- `anick.log` --- an independent combinatorial confirmation: a
  change of basis g = [[2,0,1,1], [1,0,2,0], [1,2,0,0], [0,0,1,0]] under which the four cubic initial forms
  have the combinatorially free leading multidegrees (3,3,1), (3,3,0), (3,2,2), (3,2,0), so
  Anick's criterion gives strong freeness without the Theta
  machinery.
- `mild.log` --- the Hilbert-series probe: dimensions match
  1/(1 - 4z + 4z^3) exactly through degree 6.
- `gb.log` --- a truncated Groebner basis to degree 10: 6 elements, of degrees 3, 3, 3, 3, 4, 4.
- `pre.gp`, `orakel-capped-red.gp`, `kriterium.gp` --- the scripts
  that produced the data; they run in unmodified PARI/GP.

The logs are kept verbatim as produced, including their German
labels; the verdict line "TRANSVERSALES ELEMENT GEFUNDEN => mild"
reads "transverse element found, hence mild".

## Rerunning

```sh
cd records/p3/rank-four-12144979499
DISC=-12144979499 gp -q pre.gp
POL="$(cat pol.txt)" P=3 gp -q orakel-capped-red.gp
gp -q kriterium.gp
```

The oracle writes `D-matrizen.txt` into the working directory and
needs a few minutes and a few gigabytes of PARI stack; the
criterion step reads the file from there.
