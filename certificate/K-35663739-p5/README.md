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
certified here (see `examples/p5/transverse-rank-one/`).  The
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
`examples/p5/batch-block0-01/D-35663739/result.gp` exactly; only the
certificate was kept.

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-35663739-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-35663739/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
certificate entry by entry in the notation of the paper; it is
generated deterministically from `certificate.gp` by

```sh
CERT_DIR=certificate/K-35663739-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-35663739-p5 && pdflatex certificate-guide.tex
```
