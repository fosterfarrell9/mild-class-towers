# CM class-field constructor

A deterministic alternative oracle for the unramified cyclic quintic
class fields L_x of an imaginary quadratic field: instead of the
Kummer-theoretic route through the class group and units of
K(zeta_5), the fields are constructed by complex multiplication --
double eta quotients with Schertz corrections evaluated over the
reduced forms, period sums over the five cosets of ker(x), exact
recognition, and only afterwards `polredbest`.  The construction is a
search oracle in the sense of the verification architecture: nothing
downstream trusts it, since every certificate entry passes the exact
audits and the standalone verifier.

Build (adjust the PARI prefix):

```sh
make PARI=/usr/local
```

Programs:

- `cm_construct`: enumerates the form classes, evaluates the class
  invariants at the stated precision, forms the coset period sums,
  and emits the six relative quintic polynomials with character
  metadata, the invariant, the height bound, and all rounding
  distances.
- `audit_cm_fields`: exact degree, discriminant, and Artin-kernel
  audit of a constructed field file; optionally an exact comparison
  against a stored certificate.
- `cm_character_driver`: the pipeline character driver linked with a
  wrapped `bnrclassfield` that feeds the CM fields into the otherwise
  unchanged pipeline (`--wrap=bnrclassfield`).  A builder switch is
  the intended durable interface; the wrap demonstrates
  compatibility.

Example (field of discriminant -106660295):

```sh
./cm_construct 's^2-s+26665074' cm-fields.gp 768
./audit_cm_fields cm-fields.gp ../certificates/p5/K-106660295-p5/certificate.gp
```

The certificate of the seventeenth census field, D = -109909943, was
produced through this oracle and is verified by the unchanged
standalone verifier under `certificates/p5/K-109909943-p5/`.
