# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-126740891})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 - s + 31685223 = 0
disc(K) = -126740891
Cl(K) = [195, 5, 5]
```

For this field the exhaustive GL_3(F_5) search found no Anick witness
(`MILD=UNKNOWN` in the result record); mildness is proved instead by
the transverse rank-one criterion.  The norm-degeneracy scheme has no
rational point and consists of a single reduced closed point of degree
6, whose certificate over \(\mathbf F_5[t]/(t^6+t^5+t^4+1)\) is derived
from the six matrices certified here (see
`examples/p5/transverse-rank-one/`).  The certificate pins down the
arithmetic behind the six secondary-norm matrices; everything
downstream is fast finite algebra reproducible from them.

The file contains the 18 main entries (six characters at three columns
each), exported by the per-character parallel drivers (five characters
completed with a 4 GiB stack ceiling; the character `a+c` required a
single re-run with a 16 GiB ceiling, peaking at 5.1 GiB):

```sh
cd parallelization && make PARI=/path/to/pari-prefix
python3 run_parallel.py --polynomial 's^2-s+31685223' \
  --workdir work/D-126740891 --limit exhaustive
```

Verify with the shared verifier in `certificate/`, cross-checking
against the committed result record:

```sh
cd certificate && make PARI=/path/to/pari-prefix
./verify_certificate K-126740891-p5/certificate.gp \
  ../examples/p5/further/D-126740891/result.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, `RESULT_RECORD_MATCH=PASS`, and
`CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`; from the repository root

```sh
CERT_DIR=certificate/K-126740891-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-126740891-p5 && pdflatex certificate-guide.tex
```

The repository ships the built guide only for the principal example
`certificate/K-2800905-p5`; for every other field it is these two
commands away.  The guide is exposition only; the verification chain
is implemented solely by `verify_certificate.c`.
