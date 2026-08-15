# p = 7 certificates

Arithmetic certificates for the three imaginary quadratic fields of
7-class rank three with |D_K| < 2^30, all treated in the paper:

```text
K-501510767-p7   K = Q(sqrt(-501510767)),  Cl(K) = [378, 7, 7]
K-648153647-p7   K = Q(sqrt(-648153647)),  Cl(K) = [294, 7, 7]
K-931506071-p7   K = Q(sqrt(-931506071)),  Cl(K) = [840, 7, 7]
```

All three discriminants are congruent to 1 modulo 4, so directory
names, radicands, and absolute discriminants coincide.

Each directory contains:

* `certificate.gp` — the stored arithmetic data: the six characters
  (the three coordinate characters and their three pairwise sums),
  with class fields, normalized automorphisms, norm witnesses, and
  norm classes.
* `matrices.tsv` — the secondary-norm matrices recorded by the
  search, the independent baseline the verifier's reconstruction is
  compared against.
* `transverse.txt` — the transversality protocol recomputed
  deterministically from the matrices: closed points of the
  norm-degeneracy scheme of residue degree at most three, ranks,
  reduced-isolatedness, Anick pivots, and the verdict, with a
  Singular second opinion on the point degrees.  For the first field
  three transverse points exist (one rational); for the second the
  single point of residue degree at most three has degree two and is
  transverse; for the third it is rational and transverse.

The norm witnesses of the first two fields were found through the
degree-14 class-group machinery of PARI; for the third, where a
single such computation had become impracticable, they were obtained
from the D_7-norm relations of Biasse--Fieker--Hofmann--Page
(J. London Math. Soc. 105 (2022)), computed in Hecke.  The
certificate format and the verification are identical in all three
cases.

Verify with the odd-p verifier (the prime is read from the
certificate):

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p7/K-501510767-p7/certificate.gp
```

and likewise for the other two directories.  Expected: all 18
entries verified, the shuffle identities, the six reconstructed
matrices, and `CERTIFICATE VERIFIED`.  The certificate embeds no
expected tensor (the run reports `EMBEDDED_EXPECTED_TENSOR=ABSENT`);
`matrices.tsv` holds the matrices as the search recorded them, and
comparing them with the six matrices the verifier prints is the
independent cross-check.
