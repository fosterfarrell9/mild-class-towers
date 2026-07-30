# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-123482119})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 30870530 = 0
disc(K) = -123482119
Cl(K) = [270, 5, 5]
```

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `bba, bcc, bca`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the per-character parallel drivers:

```sh
cd parallelization && make PARI=/path/to/pari-prefix
python3 run_parallel.py --polynomial 's^2-s+30870530' \
  --workdir work/D-123482119 --limit exhaustive
```

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-123482119-p5/certificate.gp \
  ../examples/p5/further/D-123482119/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
certificate entry by entry in the notation of the paper; regenerate by

```sh
CERT_DIR=certificate/K-123482119-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-123482119-p5 && pdflatex certificate-guide.tex
```
