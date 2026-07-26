# Arithmetic certificates

One subdirectory per field, named `K-<|D'|>-p5` after the radicand of the
defining polynomial.  Each subdirectory contains the data only: the
textual certificate `certificate.gp`, a README with the field data and
the exact generation command, and (for the principal example) a
human-readable guide.

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

The verifier is field-generic: it accepts 18-entry certificates (the six
main characters, as exported by the pipeline) and the principal
example's 27-entry certificate with the doubled characters.  The
verification chain and the certificate schema are documented in
`K-2800905-p5/README.md`; the required PARI 2.17.4 patch is documented
in `doc/pari-2.17.4-patch.md` at the repository root.
