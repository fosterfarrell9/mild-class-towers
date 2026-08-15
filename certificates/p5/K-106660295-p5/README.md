# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-106660295})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 26665074 = 0
disc(K) = -106660295
Cl(K) = [110, 10, 10]
```

This field has the largest class number among the twenty smallest
absolute discriminants of 5-class rank three, yet it was one of the
cheapest to compute: eleven minutes of wall time, 0.7 GiB per character
process.  The cost is governed neither by h_K nor by the class groups of
the L_x; `records/p5/README.md` records the measurement that settles
this.

The exhaustive GL_3(F_5) search found an Anick witness with leading
words `bba, bcc, bca`, so the tower group is mild (`MILD=PROVED` in the
result record).  The certificate pins down the arithmetic behind the
six secondary-norm matrices; everything downstream is fast finite
algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the per-character parallel drivers:

```sh
cd parallelization && make PARI=/path/to/pari-prefix
python3 run_parallel.py --polynomial 's^2-s+26665074' \
  --workdir work/D-106660295 --limit exhaustive
```

Verify with the shared verifier in `certificates/p5/`, cross-checking
against the committed result record:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p5/K-106660295-p5/certificate.gp \
  ../records/p5/further/D-106660295/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificates/p5/K-106660295-p5 gp -qf tools/certificate_guide.gp
cd certificates/p5/K-106660295-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificates/p5/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
