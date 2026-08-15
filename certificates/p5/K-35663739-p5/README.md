# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-35663739})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 8915935 = 0
disc(K) = -35663739
Cl(K) = [30, 10, 5]
```

For this field the exhaustive GL_3(F_5) search found no Anick witness
(`MILD=UNKNOWN` in the result record); mildness is proved instead by
the transverse rank-one criterion, whose degree-6 certificate over
\(\mathbf F_5[t]/(t^6+t^5+t^4+1)\) is derived from the six matrices
certified here (see `records/p5/transverse-rank-one/`).  The
certificate pins down the arithmetic behind the six secondary-norm
matrices; everything downstream is fast finite algebra reproducible
from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the batch runner's certificate re-run mode on the
already computed field:

```sh
python3 tools/run_mildness_batch.py --resume \
  --certificate-dir certificate --certificates --certificates-only
```

The fresh run reproduced the verified entries of the committed
`records/p5/batch-block0-01/D-35663739/result.gp` exactly; only the
certificate was kept.

Verify with the shared verifier in `certificates/p5/`, cross-checking
against the committed result record:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p5/K-35663739-p5/certificate.gp \
  ../records/p5/batch-block0-01/D-35663739/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificates/p5/K-35663739-p5 gp -qf tools/certificate_guide.gp
cd certificates/p5/K-35663739-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificates/p5/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
