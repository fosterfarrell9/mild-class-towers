# p = 7 certificates

Arithmetic certificates for the two imaginary quadratic fields of
7-class rank three treated in the paper:

```text
K-501510767-p7   K = Q(sqrt(-501510767)),  Cl(K) = [378, 7, 7]
K-648153647-p7   K = Q(sqrt(-648153647)),  Cl(K) = [294, 7, 7]
```

Both discriminants are congruent to 1 modulo 4, so directory names,
radicands, and absolute discriminants coincide.

Each directory contains:

* `certificate.gp` — the stored arithmetic data: six characters and
  three doubled characters, with class fields, normalized
  automorphisms, auxiliary pairs, and norm classes.
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
  transverse.

Verify with the odd-p verifier (the prime is read from the
certificate):

```sh
make -C examples/p3/verifier PARI=/path/to/pari-prefix
examples/p3/verifier/verify_certificate examples/p7/K-501510767-p7/certificate.gp
examples/p7/K-648153647-p7/certificate.gp likewise
```

Expected: all 18 entries verified, the doubled-character identities
D_{2x} = 4 D_x, the shuffle identities, agreement of the
reconstructed tensor with `matrices.tsv`, and `CERTIFICATE
VERIFIED`.
