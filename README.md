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
- `records/p3/rank-four-653329427/` records the single field of the range
  with 3-class rank four, reported as undecided in Section 5.5 of the
  paper: the ten secondary norm operators, the Bockstein matrix of rank
  four and the empty cone, with the pure-GP scripts and an independent
  second run.
- Three fields of 3-class rank four beyond the range of the paper are
  reported mild, with arithmetic certificates extending the paper's
  verification layer to rank four; see the next section.

## 🎉 Beyond the paper: three mild fields of 3-class rank four (announced 26 August 2026; certified 27 August 2026)

A sweep of the fields of 3-class rank four with |D_K| < 20294967296,
run with the arithmetic of Section 4 of the paper adapted to rank
four, found three fields whose 3-class tower group satisfies the
transversality criterion --- each at a single rational point of its
Bockstein cone, over F_3 itself:

| D_K | class group | 3-class group | rk D_x | rk Theta_x |
| --- | --- | --- | --- | --- |
| -12144979499 | [1278, 3, 3, 3] | (9, 3, 3, 3) | 2 = d-2 | 2 (surjective) |
| -18191474648 | [1152, 6, 3, 3] | (9, 3, 3, 3) | 2 = d-2 | 2 (surjective) |
| -18561189299 | [1674, 3, 3, 3] | (9, 3, 3, 3) | 2 = d-2 | 2 (surjective) |

By the criterion (Theorem B of the paper, stated for any rank
d >= 3) each of these groups is mild and has cohomological
dimension 2.  To our knowledge, they are the first mild p-class
tower groups of p-class rank four.

Each field carries an arithmetic certificate in the sense of the
paper's verification appendix, extended to rank four: the standard
character family with explicit norm witnesses, 40 entries per
field, built by `tools/build_witness_certificate.gp` and checked by
the independent pure-GP verifier
`tools/verify_certificate_general.gp` (812 checks per field; ideal
arithmetic plus the base class group made unconditional with
bnfcertify).  The reconstructed tensors agree entry for entry with
the announced records, and the transversality verdicts on them are
exact linear algebra over F_3.  Each computation was also
reproduced from scratch with byte-identical results and confirmed
by three independent probes (Anick's criterion, the Hilbert series,
a truncated Groebner basis).

The full data --- tensors, criterion runs, crosscheck runs, probe
logs, and the scripts, which run in unmodified PARI/GP --- are in
`records/p3/rank-four-12144979499/`,
`records/p3/rank-four-18191474648/`, and
`records/p3/rank-four-18561189299/`.

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
`doc/pari-2.17.4-patch.md`.  No PARI source is redistributed here.  A single
certificate can be checked without that build, in unmodified PARI/GP, with
`tools/verify_certificate.gp`; see `REPRODUCING.md`.

## License

MIT, see `LICENSE`.  Copyright is held by Eric Ahlqvist for the original
program and by Denis Vogel for the work added since; `NOTICE` lists which
files belong to which part.
