# Certificate for D = -198040904

Imaginary quadratic field of discriminant `-198040904`, class group `[396,9,3]`,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../../../records/p3/FORMAT.md` for the
format and `../../../../records/p3/README.md` for the workflow):

```sh
make -C ../../../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/019/K-198040904-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`1bd73ed0e60e5700adab688121f1e38d4a7974c59a9606785e3be8b5b4bf7a60`
