# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-235466845})\), \(p=5\)

Certificate for the secondary-norm matrices of the field with

```text
p = 5
K = Q(s), s^2 + 235466845 = 0
disc(K) = -941867380
Cl(K) = [180,10,5]
```

(The directory is named after the radicand 235466845; the
discriminant is \(4\cdot(-235466845)=-941867380\).)

The file contains the 18 main entries (six characters at three
columns each), produced by the slice runners and re-verified by the
standalone verifier on a second machine.

Verify:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificate/K-235466845-p5/certificate.gp
```

Expected: `BASE_BNF_CERTIFIED=PASS`, all 18 entry lines with
`AC1=PASS AC2=PASS`, the six matrices, and `CERTIFICATE VERIFIED`.

A human-readable guide restating the certificate entry by entry in
the notation of the paper is generated deterministically from
`certificate.gp`:

```sh
CERT_DIR=certificate/K-235466845-p5 gp -qf tools/certificate_guide.gp
cd certificate/K-235466845-p5 && pdflatex certificate-guide.tex
```
