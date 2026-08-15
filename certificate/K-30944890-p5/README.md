# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-30944890})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 + 30944890 = 0
disc(K) = -123779560
Cl(K) = [10, 10, 10, 2]
```

(The directory is named after the radicand 30944890; the discriminant
is \(4\cdot(-30944890)=-123779560\).)

For this field the exhaustive GL_3(F_5) search found an Anick witness
with leading words `cca, cbb, cba`, so the tower group is mild
(`MILD=PROVED` in the result record).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything downstream
is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the per-character parallel drivers:

```sh
cd parallelization && make PARI=/path/to/pari-prefix
python3 run_parallel.py --polynomial 's^2+30944890' \
  --workdir work/D-123779560 --limit exhaustive
```

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificate/K-30944890-p5/certificate.gp \
  ../examples/p5/further/D-123779560/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificate/K-30944890-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-30944890-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificate/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
