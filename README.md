# mild-class-towers

Companion repository for the paper *Mild p-class tower groups of
imaginary quadratic fields*.  It is the continuation of Eric
Ahlqvist's `Massey-pari`; the authorship boundary is recorded in
`NOTICE` and visible directly in the git history.  It contains

- arithmetic certificates for all 12 956 imaginary quadratic fields of
  odd p-class rank three with |D_K| < 2^30 — 12 749 at p = 3, 204 at
  p = 5, 3 at p = 7 — and the standalone verifier (`verifier/`) that
  checks each certificate by exact arithmetic, independently of the
  search that produced it, together with arithmetic witness
  certificates for the fields of 3-class rank four through 2 * 10^10;
- the transversality, strong-freeness, and cone-criterion records
  behind the mildness results, the flag certificates of Part W7, and
  the parts W1–W7 cited by the paper;
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

## Relation to the paper

This repository accompanies the paper *Mild p-class tower groups of
imaginary quadratic fields*,
[arXiv:2608.17072](https://arxiv.org/abs/2608.17072).  The text this
repository documents is
[version 3](https://arxiv.org/abs/2608.17072v3); every part the paper
cites — the parts W1–W7 under `worksheets/`, the certificate
collections, the records and the verifiers — is in this tree.
Releases are archived on Zenodo; the initial release v1.0.0 has
[doi:10.5281/zenodo.21982734](https://doi.org/10.5281/zenodo.21982734).

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
| `records/p3/flag-certificates/` | Part W7 (rank three, p = 3): one certificate per field decided by a word count after a change of variables, with the change of variables, the coefficient field and the high terms; format note, worked examples, index |
| `records/p5/flag-certificates/` | Part W7 (p = 5): the certificate for the field D = -781922404, with its relation tensor under `records/p5/D-781922404/` |
| `records/p3/flag-certificates-rank-four/` | Part W7 (rank four): the certificate for the rank-four field D = -653329427; its witness certificate, verification log and tensor are under `records/p3/rank-four-653329427/` |
| `records/p3/rank-four-census/` | the rank-four census through 2 * 10^10 plus the next field D = -20217903567: manifest and, per field, the relation tensor, the flag certificate, the witness certificate and its verification log; the three fields decided by the transversality criterion keep their standalone records under `records/p3/rank-four-12144979499/`, `-18191474648/`, `-18561189299/` |
| `census/` | the enumeration behind the census: the complete discriminant lists per prime, with the Mosunov--Jacobson provenance |
| `worksheets/` | the parts cited by the paper, named by part: `W1-W4-worksheets.pdf`, `W5-tables-p5.pdf`, `W6-tables-p3.pdf`, and `part-W7.tsv`, the index of the fields decided by a word count after a change of variables (see `worksheets/README.md`) |
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
`doc/pari-2.17.4-patch.md`.  No PARI source is redistributed here.  A single
certificate can be checked without that build, in unmodified PARI/GP, with
`tools/verify_certificate.gp`; see `REPRODUCING.md`.

## License

MIT, see `LICENSE`.  Copyright is held by Eric Ahlqvist for the original
program and by Denis Vogel for the work added since; `NOTICE` lists which
files belong to which part.
