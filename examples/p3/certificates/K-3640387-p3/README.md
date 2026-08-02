# Certificate for D = -3640387

Imaginary quadratic field of discriminant `-3640387`, class group ``,
prime `p = 3`.  Role: `pilot`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
../../verifier/verify_certificate certificate.gp
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`73a76d685c5de83d14c97e641c65311ccb538de35490b37b6ccfd61b7643d0b1`
