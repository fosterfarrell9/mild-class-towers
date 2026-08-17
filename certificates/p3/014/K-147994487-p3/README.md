# Certificate for D = -147994487

Imaginary quadratic field of discriminant `-147994487`, class group `[432,9,3]`,
prime `p = 3`.  Role: `theorem_7_1`.

Verify with the standalone verifier (see `../../../../records/p3/FORMAT.md` for the
format and `../../../../records/p3/README.md` for the workflow):

```sh
make -C ../../../../verifier PARI=/usr/local
verifier/verify_certificate certificates/p3/014/K-147994487-p3/certificate.gp   # from the repository root
```

Expected markers: eighteen entries each ending in `AC1=PASS` and
`AC2=PASS`, then `SHUFFLE_IDENTITIES=PASS`,
`EXPECTED_TENSOR_MATCH=PASS`, and the final line
`CERTIFICATE VERIFIED`.

SHA-256 of `certificate.gp`:
`52300a54f72c0bb34f6e36d0871d6d0a5859d25fa5615077063ece8357b0517e`
