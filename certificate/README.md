# Arithmetic certificates

One subdirectory per field, named `K-<|D'|>-p5` after the radicand of the
defining polynomial.  Each subdirectory contains the data only: the
textual certificate `certificate.gp`, a README with the field data and
the exact generation command, and a human-readable guide
`certificate-guide.pdf` generated from the certificate by
`tools/certificate_guide.gp`.

The verifier is shared and lives here, next to this file.  Build it once
and point it at any certificate:

```sh
cd certificate
make PARI=/path/to/pari-prefix
./verify_certificate K-2800905-p5/certificate.gp
./verify_certificate K-51213139-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-51213139/result.gp
```

The optional second argument names a committed `result.gp` record; the
verifier then cross-checks the reconstructed matrices against its
`secondary_norm_samples` and reports `RESULT_RECORD_MATCH=PASS`.

## Integral bases, and why the format has a version

Elements and ideals are stored as coordinates with respect to an integral
basis, and the basis PARI returns is LLL-reduced -- hence not canonical.
Another machine's `nfinit` may reduce differently, and the same coordinate
vector then denotes a different algebraic number: the stored automorphism
stops fixing `K`, the stored ideals become other ideals, and verification
fails on data that is entirely correct.  This is not hypothetical.  When
the twenty-one certificates were first checked on a second machine, eleven
of them failed or ran out of time.

Format 2 therefore records the basis it was written against: in the header
for `K`, and as the last field of every entry for the class field of that
entry.  The verifier expresses the stored basis in its own, requires the
resulting matrix to be integral with determinant `±1` -- so that the two
bases span the same ring of integers, and a certificate cannot substitute
a different order -- and converts the coordinates before checking anything.
Where the bases agree the matrix is the identity.

`upgrade_certificate` converts a format 1 certificate, and must be run on
a machine where that certificate verifies: only there does `nfinit` return
the basis its coordinates were written against.  It recomputes no
arithmetic.

```sh
./upgrade_certificate old/certificate.gp new-certificate.gp
```

The verifier still reads format 1, with a warning on stderr, since a
failure then carries no information about the arithmetic.

## Scope

The verifier is field-generic: it accepts 18-entry certificates (the six
main characters, as exported by the pipeline) and the principal
example's 27-entry certificate with the doubled characters.  The
verification chain and the certificate schema are documented in
`K-2800905-p5/README.md`; the required PARI 2.17.4 patch is documented
in `doc/pari-2.17.4-patch.md` at the repository root.
