# Massey-pari

A C program using the [PARI](http://pari.math.u-bordeaux.fr/) library that
computes Massey products in the étale cohomology of the ring of integers of
a number field, together with the arithmetic that goes with them.  PARI
makes it fast for number fields of low degree.

The original program was written by Eric Ahlqvist for the computations in
Ahlqvist--Carlson, [*Massey products in the étale cohomology of number
fields*](https://www.degruyter.com/document/doi/10.1515/crelle-2025-0006/html).
This repository continues that work: it adds prescribed-character secondary
norm operators, the reconstruction of the complete cubic relation space,
strong-freeness and mildness tests, and a certificate infrastructure by
which the arithmetic can be checked independently of the search that
produced it.  `NOTICE` records the authorship boundary, which is also
visible directly in the git history.

## Quick start

```sh
tools/build-patched-pari.sh            # patched PARI 2.17.4 into ~/.local
make PARI="$HOME/.local"               # the main program
make -C verifier PARI="$HOME/.local"
verifier/verify_certificate certificate/K-2800905-p5/certificate.gp
```

The last command reverifies, by exact arithmetic and without repeating any
search, the secondary norm operators of the principal example.  The original
program is run as

```sh
./build/massey p "pol(s)"
```

with `pol(s)` a defining polynomial of the number field K in the variable
`s`, and `p` a prime dividing the class number.

## Layout

| Directory | Contents |
|---|---|
| `src/`, `headers/` | the C program: field construction, Artin symbols, secondary norm operators, Massey tensor, cubic relation space |
| `verifier/` | the standalone certificate verifier shared by all three primes, with its Makefile and rejection tests |
| `certificate/` | one directory per certified p = 5 field, with a human-readable guide per certificate |
| `examples/p5/` | computed results and run logs per field, and the transverse rank-one certificates |
| `examples/p3/` | the p = 3 collection: certificates and source tensors in discriminant buckets, result records, cone-criterion reports, and drivers |
| `examples/p7/` | the three p = 7 fields: certificates, matrices, and transversality protocols |
| `worksheets/` | the worksheet parts W1--W6 cited by the paper |
| `tools/` | batch runner, Gröbner and Singular drivers, transverse rank-one certification, class-group table readers, PARI build script |
| `parallelization/` | per-character parallel drivers for expensive fields |
| `tests/` | regression tests |
| `doc/` | the required one-file PARI patch |

`REPRODUCING.md` documents every computation and how to rerun it.

## Certificates

For each computed field the repository stores a textual, PARI-readable
*arithmetic certificate*: the class field, its normalized automorphism, the
Ahlqvist--Carlson data, and the resulting norm classes.  The standalone
verifier recomputes and certifies the base field, identifies the class field
and its generator from the stored data, checks the Ahlqvist--Carlson
identities by exact ideal arithmetic, and recovers each norm class — without
rerunning the search, and without computing any class or unit group of the
relative fields.  A successful run ends with `CERTIFICATE VERIFIED`.

Building the verifier requires the same patched PARI; the patch makes three
internal routines non-static and is documented in
`doc/pari-2.17.4-patch.md`.  No PARI source is redistributed here.

## License

MIT, see `LICENSE`.  Copyright is held by Eric Ahlqvist for the original
program and by Denis Vogel for the work added since; `NOTICE` lists which
files belong to which part.
