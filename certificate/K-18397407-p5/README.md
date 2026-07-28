# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-18397407})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 4599352 = 0
disc(K) = -18397407
Cl(K) = [40, 10, 5]
```

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `bba, bcc, bca`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by re-running the audited pipeline on the already
computed field:

```sh
MASSEY_CERTIFICATE_EXPORT="$PWD/certificate/K-18397407-p5/certificate.gp" \
./build/massey --example-result /tmp/result.gp \
  --strong-search-limit exhaustive 5 's^2-s+4599352'
```

The fresh run reproduced the verified entries of the committed
`examples/p5/K-18397407/result.gp` exactly; only the certificate was
kept.

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-18397407-p5/certificate.gp \
  ../examples/p5/K-18397407/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
certificate entry by entry in the notation of the paper; it is
generated deterministically from `certificate.gp` by

```sh
CERT_DIR=certificate/K-18397407-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-18397407-p5 && pdflatex certificate-guide.tex
```
