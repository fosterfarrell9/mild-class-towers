# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-13579778})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 + 13579778 = 0
disc(K) = -54319112
Cl(K) = [60, 10, 5]
```

(The directory is named after the radicand 13579778; the discriminant
is \(4\cdot(-13579778)=-54319112\).)

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `cca, cbb, cba`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the batch runner's certificate re-run mode on the
already computed field:

```sh
python3 tools/run_mildness_batch.py --resume \
  --certificate-dir certificate --certificates --certificates-only
```

The fresh run reproduced the verified entries of the committed
`examples/p5/batch-block0-01/D-54319112/result.gp` exactly; only the
certificate was kept.

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificate/K-13579778-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-54319112/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificate/K-13579778-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-13579778-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificate/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
