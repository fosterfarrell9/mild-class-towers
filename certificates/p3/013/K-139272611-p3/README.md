# Certificate for D = -139272611

Imaginary quadratic field of discriminant `-139272611`, class group ``,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/013/K-139272611-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`18c325e74537172e52e3d0610b47d7b779739315eb1b2c602be7c8aeb300e7ab`
