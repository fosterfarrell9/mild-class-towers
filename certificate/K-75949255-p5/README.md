# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-75949255})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 18987314 = 0
disc(K) = -75949255
Cl(K) = [40, 20, 5]
```

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `ccb, caa, cab`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
(the cubic matrix T, the witness check) is fast finite algebra
reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported inline by the batch run

```sh
python3 tools/run_mildness_batch.py --resume \
  --certificate-dir certificate --limit 1
```

(equivalently, by the audited pipeline with
`MASSEY_CERTIFICATE_EXPORT` set; see `REPRODUCING.md` at the repository
root).

Verify with the shared verifier in `certificate/`, optionally
cross-checking against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-75949255-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-75949255/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
certificate entry by entry in the notation of the paper; it is
generated deterministically from `certificate.gp` by

```sh
CERT_DIR=certificate/K-75949255-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-75949255-p5 && pdflatex certificate-guide.tex
```
