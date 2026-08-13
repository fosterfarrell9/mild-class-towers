# Arithmetic certificates

One subdirectory per field, named `K-<|D'|>-p5` after the radicand of the
defining polynomial.  Each subdirectory contains the data only: the
textual certificate `certificate.gp` and a README with the field data and
the exact generation command.

A human-readable guide `certificate-guide.pdf` restating a certificate
entry by entry in the notation of the paper is generated from
`certificate.gp` by `tools/certificate_guide.gp`; from the repository
root

```sh
CERT_DIR=certificate/K-<n>-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-<n>-p5 && pdflatex certificate-guide.tex
```

The built guide ships only for the principal example `K-2800905-p5`,
since it is derived from data that ships in full.  The guide is
exposition; the verification chain is implemented solely by
`verify_certificate.c`.

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

`make check` first verifies `K-2800905-p5/certificate.gp` unchanged and
then runs `test_rejections.py`, which generates temporary copies of
that certificate and requires the verifier to reject every one of them
with the expected message: five copies alter the arithmetic content ---
the character vector, the normalized automorphism, the multiplicity of
an entry, a norm-class vector, and the absolute field model, each
separately --- and two alter the container (an unsupported format
version, a missing entry).  The outcomes are recorded in
`rejection-tests.json`.

## Integral bases

Elements and ideals are stored as coordinates with respect to an integral
basis of the field in question.  PARI's integral basis is LLL-reduced and
therefore not canonical: another `nfinit` may return a different one, and
the same coordinate vector would then denote a different algebraic number.

Each certificate therefore records the bases it refers to: in the header
for `K`, and as the last field of every entry for the class field of that
entry.  The verifier expresses the stored basis in its own, requires the
resulting matrix to be integral with determinant `±1` -- so that both bases
span the same ring of integers, and a certificate cannot substitute a
different order -- and converts the coordinates before checking anything.
Where the two bases agree, that matrix is the identity.

To check the whole collection at once, from the repository root:

```sh
python3 tools/verify_all_certificates.py
```

It prints one line per field and exits nonzero unless every certificate
verified.  Running it on a machine that produced none of them is what
tests the arrangement above.

## Scope

The verifier is field-generic: it accepts 18-entry certificates (the six
main characters, as exported by the pipeline) and the principal
example's 27-entry certificate with the doubled characters.  The
verification chain and the certificate schema are documented in
`K-2800905-p5/README.md`; the required PARI 2.17.4 patch is documented
in `doc/pari-2.17.4-patch.md` at the repository root.
