# A mild field of 3-class rank four (certified)

The imaginary quadratic field with D_K = -18191474648 has class group of
invariants [1152, 6, 3, 3], class number 62208, and elementary quotient of
3-class group (9, 3, 3, 3).  Its Bockstein cone is a single
rational point, and at that point both parts of the transversality
criterion hold:

```
rk D_x = 2 = d-2            (rank condition)
Theta matrix = [1, 2; 1, 0; 1, 0],  rank 2  (surjective)
```

The point is a transverse element over F_3 itself, so by the
transversality criterion of the paper the 3-class tower group of
the field is mild, of cohomological dimension 2; see the status below.  This is one of three such fields found on
26 August 2026; see the repository README.

## Status

Certified 27 August 2026.  `certificate.gp` is an arithmetic
certificate in the sense of the paper's verification appendix,
extended to rank four: the standard character family x_i, x_i+x_j,
with 40 entries carrying explicit norm witnesses (t, I').
`verification.log` records its check by
`tools/verify_certificate_general.gp`, an independent pure-GP
verifier (812 checks; ideal arithmetic plus the base class group
made unconditional with bnfcertify -- the only certified class
group the chain needs).  The reconstructed tensor agrees entry for
entry with `D-matrizen.txt`.  Before certification the runs had
been reproduced from scratch with byte-identical results and
confirmed by the three independent probes below.

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
  change of basis g = [[2,2,2,1], [2,0,0,0], [0,2,1,0], [1,2,0,0]] under which the four cubic initial forms
  have the combinatorially free leading multidegrees (3,3,1), (3,3,0), (3,2,2), (3,2,0), so
  Anick's criterion gives strong freeness without the Theta
  machinery.
- `mild.log` --- the Hilbert-series probe: dimensions match
  1/(1 - 4z + 4z^3) exactly through degree 6.
- `gb.log` --- a truncated Groebner basis to degree 10: 6 elements, of degrees 3, 3, 3, 3, 4, 4.
- `certificate.gp` --- the arithmetic certificate: for each of the
  ten standard characters and each of the four torsion columns the
  field model, sigma, the pair (a', J), the norm witness (t, I'), a
  sign-fixing prime, and the norm class.  Built by
  `tools/build_witness_certificate.gp`.
- `build.log`, `verification.log` --- the certificate build run and
  the full verifier report.
- `pre.gp`, `orakel-capped-red.gp`, `kriterium.gp` --- the scripts
  that produced the search data; they run in unmodified PARI/GP.

The logs are kept verbatim as produced, including their German
labels; the verdict line "TRANSVERSALES ELEMENT GEFUNDEN => mild"
reads "transverse element found, hence mild".

## Rerunning

```sh
cd records/p3/rank-four-18191474648
DISC=-18191474648 gp -q pre.gp
POL="$(cat pol.txt)" P=3 gp -q orakel-capped-red.gp
gp -q kriterium.gp
```

The oracle writes `D-matrizen.txt` into the working directory and
needs a few minutes and a few gigabytes of PARI stack; the
criterion step reads the file from there.  To recheck the
certificate:

```sh
CERT_DIR=. gp -q ../../../tools/verify_certificate_general.gp
```
