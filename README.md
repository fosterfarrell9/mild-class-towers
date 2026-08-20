# mild-class-towers

Companion repository for the paper *Mild p-class tower groups of
imaginary quadratic fields*.  It is the continuation of Eric
Ahlqvist's `Massey-pari`; the authorship boundary is recorded in
`NOTICE` and visible directly in the git history.  It contains

- arithmetic certificates for all 12 956 imaginary quadratic fields of
  odd p-class rank three with |D_K| < 2^30 — 12 749 at p = 3, 204 at
  p = 5, 3 at p = 7 — and the standalone verifier (`verifier/`) that
  checks each certificate by exact arithmetic, independently of the
  search that produced it;
- the transversality, strong-freeness, and cone-criterion records
  behind the mildness results, and the worksheet parts W1–W6 cited in
  the paper;
- the C/[PARI](http://pari.math.u-bordeaux.fr/) search program that
  computes Massey products in the étale cohomology of rings of
  integers of number fields, with prescribed-character secondary norm
  operators, the reconstruction of the complete cubic relation space,
  and the strong-freeness and mildness tests.

The original program was written by Eric Ahlqvist for the computations
in Ahlqvist--Carlson, [*Massey products in the étale cohomology of
number
fields*](https://www.degruyter.com/document/doi/10.1515/crelle-2025-0006/html);
the tag `ahlqvist-final` marks his last commit.

## Relation to the published paper

The paper cites version **v1.0.0**, archived at
[doi:10.5281/zenodo.21982734](https://doi.org/10.5281/zenodo.21982734).  Every
number and every claim there refers to that tag.  The certificates, the records
and the results are unchanged since.

The current state differs from that tag in the following respects:

- A single certificate can be checked in unmodified PARI/GP, with
  `tools/verify_certificate.gp`.  The patched PARI is needed only for the
  sweep over all certificates.

## Quick start

```sh
tools/build-patched-pari.sh            # patched PARI 2.17.4 into ~/.local
make PARI="$HOME/.local"               # the main program
make -C verifier PARI="$HOME/.local"
verifier/verify_certificate certificates/p5/K-2800905-p5/certificate.gp
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
| `certificates/` | the arithmetic certificates: one collection per prime (`p3/` in discriminant buckets, `p5/`, `p7/`), one directory per field |
| `records/` | what the searches produced and what follows from the verified matrices: result records and run logs, source tensors, verification and strong-freeness records, cone-criterion reports, transversality certificates |
| `census/` | the enumeration behind the census: the complete discriminant lists per prime, with the Mosunov--Jacobson provenance |
| `worksheets/` | the worksheet parts W1--W6 cited by the paper |
| `tools/` | the drivers and generators: verification harnesses, batch runner, Gröbner and Singular engines, transverse rank-one certification, class-group table readers, PARI build script |
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
