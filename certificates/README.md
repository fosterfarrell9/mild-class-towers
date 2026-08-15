# Arithmetic certificates

One collection per prime, one directory per field — 12 956 fields in
total, the complete census of imaginary quadratic fields of odd
p-class rank three with |D_K| < 2^30:

- `p3/<bucket>/K-<|D|>-p3/` — the 12 749 fields at p = 3, sharded
  into buckets `<bucket>` = floor(|D|/10^7) as three digits
  (`000`..`107`);
- `p5/K-<m>-p5/` — the 204 fields at p = 5, named after the radicand
  `m` of the defining polynomial `s^2 + m` (so |D_K| = m for
  m = 3 mod 4 and |D_K| = 4m otherwise);
- `p7/K-<|D|>-p7/` — the three fields at p = 7 (radicand and |D_K|
  coincide).

Every field directory contains the textual data file
`certificate.gp`.  Some p = 3 directories carry a `hints.gp` sidecar
of proven prime factors (see `records/p3/README.md`); the p = 5
directories carry a README with the field data, and the p = 7
directories additionally hold the search matrices and transversality
protocols (see `p7/README.md`).

The shared verifier lives in the top-level directory `verifier/`;
from the repository root:

```sh
make -C verifier PARI=/path/to/pari-prefix
verifier/verify_certificate certificates/p5/K-2800905-p5/certificate.gp
```

The certificate format and the verification chain are documented in
Appendix C of the paper, in `records/p3/FORMAT.md`, and in
`p5/README.md`.
