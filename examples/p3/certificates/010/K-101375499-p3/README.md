# Certificate for D = -101375499

Imaginary quadratic field of discriminant `-101375499`, class group ``,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
verifier/verify_certificate examples/p3/certificates/010/K-101375499-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`bbb59a6a863a67743a2ac1e387d1b563e16c7ad8eb50929d6d2cecc6325f3477`
