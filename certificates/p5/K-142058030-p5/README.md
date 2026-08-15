# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-142058030})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 + 142058030 = 0
disc(K) = -568232120
Cl(K) = [50,10,10,2]
```

(The directory is named after the radicand 142058030; the
discriminant is \(4\cdot(-142058030)=-568232120\).)

The file contains the 18 main entries (six characters at three
columns each), produced by the slice runners and re-verified by the
standalone verifier on a second machine.

Verify:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p5/K-142058030-p5/certificate.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, and `CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`:

```sh
CERT_DIR=certificates/p5/K-142058030-p5 gp -qf tools/certificate_guide.gp
cd certificates/p5/K-142058030-p5 && pdflatex certificate-guide.tex
```
