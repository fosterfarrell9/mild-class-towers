# Certificate for D = -3640387

Imaginary quadratic field of discriminant `-3640387`, class group `[18,3,3]`,
prime `p = 3`.  Role: `pilot`.

Verify with the standalone verifier (see `../../../../records/p3/FORMAT.md` for the
format and `../../../../records/p3/README.md` for the workflow):

```sh
make -C ../../../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/000/K-3640387-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`73a76d685c5de83d14c97e641c65311ccb538de35490b37b6ccfd61b7643d0b1`
