# Certificate for D = -228404408

Imaginary quadratic field of discriminant `-228404408`, class group ``,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../FORMAT.md` for the
format and `../../README.md` for the workflow):

```sh
make -C ../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/022/K-228404408-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`2b795f2de59ddf22f3c903d505c4d7eb3b8885630c7c587e316625e0c884570a`
