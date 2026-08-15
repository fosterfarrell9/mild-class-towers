# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-77778287})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 19444572 = 0
disc(K) = -77778287
Cl(K) = [120, 10, 5]
```

This field was the most expensive of the twenty smallest absolute
discriminants of 5-class rank three, and for a while it looked
intractable: three runs died at ceilings of 4, 8 and 12 GiB, the last of
them minutes before it would have finished.  The obstruction was memory
alone.  Every character needs 16.75 GiB, and once that was granted each
one completed -- character `a` in 9 h 14 min here, the remaining five in
about 1 h 42 min each on a faster machine.

The six peaks are 17564236, 17563508, 17563328, 17563264, 17563272 and
17563800 kbytes: a spread of 0.006 %.  The memory is therefore consumed
by a computation that does not depend on the character at all, which is
why splitting a field across six processes buys wall time but no
headroom.

The exhaustive GL_3(F_5) search found an Anick witness with leading
words `bba, bcc, bca`, so the tower group is mild (`MILD PASS` in the
result record).  The certificate pins down the arithmetic behind the six
secondary-norm matrices; everything downstream is fast finite algebra
reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the per-character parallel drivers:

```sh
cd parallelization && make PARI=/path/to/pari-prefix
MASSEY_PARISTACK_MAX=25769803776 python3 run_parallel.py \
  --polynomial 's^2-s+19444572' --workdir work/D-77778287 \
  --limit exhaustive
```

The ceiling matters: with the default of 8 GiB the run dies partway
through.  On a machine that cannot hold six such processes at once, use
`--characters` to compute them in waves and `--merge-only` to assemble
the certificate afterwards.

Verify with the shared verifier in `certificates/p5/`, cross-checking
against the committed result record:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p5/K-77778287-p5/certificate.gp \
  ../records/p5/further/D-77778287/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.  The contrast with the search is the point of the
certificate: nine seconds on the machine that took nine hours per
character, because verification repeats no search and computes no class
or unit group of the relative fields.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificates/p5/K-77778287-p5 gp -qf tools/certificate_guide.gp
cd certificates/p5/K-77778287-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificates/p5/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
