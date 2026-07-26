# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-51213139})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 12803285 = 0
disc(K) = -51213139
Cl(K) = [75, 5, 5]
```

one of the two fields of the paper's appendix whose mildness status is
open.  The certificate pins down the arithmetic behind the six
secondary-norm matrices; everything downstream (the cubic matrix T, the
exhaustive GL_3(F_5) search, the Groebner/Hilbert evidence in
`examples/p5/batch-block0-01/D-51213139/strong-freeness/`) is fast
finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the audited example pipeline; see `REPRODUCING.md`
at the repository root.  It was generated with

```sh
MASSEY_CERTIFICATE_EXPORT=certificate/K-51213139-p5/certificate.gp \
./build/massey --example-result /tmp/result.gp \
  --strong-search-limit exhaustive 5 's^2-s+12803285'
```

Verify with the shared verifier in `certificate/`, optionally
cross-checking against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-51213139-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-51213139/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
certificate entry by entry in the notation of the paper; it is
generated deterministically from `certificate.gp` by

```sh
CERT_DIR=certificate/K-51213139-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-51213139-p5 && pdflatex certificate-guide.tex
```
