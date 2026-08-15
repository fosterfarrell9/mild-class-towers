# Certificate for D = -138230347

Imaginary quadratic field of discriminant `-138230347`, class group ``,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/013/K-138230347-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`928ef0825747aeec56a6a38fe6cc897b87a86e835894948084d0718b16819324`
