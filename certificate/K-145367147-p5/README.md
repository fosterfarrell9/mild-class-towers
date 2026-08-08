# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-145367147})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 36341787 = 0
disc(K) = -145367147
Cl(K) = [125, 5, 5]
```

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `aac, abb, abc`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each).  It was exported by the per-character parallel drivers in
`parallelization/` (six audited single-character processes; the six
partial certificates were merged in the canonical entry order with
identical headers):

```sh
cd parallelization && make PARI=/path/to/pari-prefix
python3 run_parallel.py --polynomial 's^2-s+36341787' \
  --workdir work/D-145367147 --limit exhaustive
```

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-145367147-p5/certificate.gp \
  ../examples/p5/batch-block0-01/D-145367147/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificate/K-145367147-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-145367147-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificate/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
