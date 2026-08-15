# Certificate for D = -53209523

Imaginary quadratic field of discriminant `-53209523`, class group ``,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
verifier/verify_certificate examples/p3/certificates/005/K-53209523-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`04b6720bb7da423691571cf2fbaa4bd4702e4ee481e2ec11f82d3b391b05d2bf`
