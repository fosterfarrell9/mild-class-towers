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

## 🎉 Beyond the paper: 11 379 more mild fields of 3-class rank three, the last field at p = 5, and the complete rank-four census (announced 1 September 2026)

The paper left the mildness of 11 859 fields of the 3-rank-3 census
undecided.  For 11 379 of them the question is now settled: each of
these fields has a mild 3-class tower group of cohomological
dimension 2.

The proofs are certificates of the kind used in the paper --- a finite
Groebner basis of the cubic relation ideal whose word counts match the
series 1/(1-3z+3z^3) (Appendix A of the paper) --- computed after an
invertible linear change of variables: for 9 975 fields over F_3
itself, for 1 404 fields over an extension field F_{3^e} with
2 <= e <= 8 (811 over F_9, 277 over F_27, 152 over F_81, 81 over
F_243, 42 over F_729, 15 over F_2187, 26 over F_6561).  A change of
variables over an extension field suffices, because the Hilbert
series is invariant under base change; `records/p3/flag-certificates/
FORMAT.md` gives the complete argument, which uses only statements of
the paper.

Every one of the 11 379 certificates passes the standalone verifier
`tools/verify_flag_certificate.py` (plain Python, no computer algebra
system: it recomputes the Groebner basis and the word counts from the
relation tensors of the census records), and every one has been
recomputed independently by a second engine --- the certificates over
F_3 by the Groebner engine of the paper (`tools/strong-freeness`), the
certificates over extension fields by Singular's Letterplace
implementation with the certificate's minimal polynomial as `minpoly`.
Three certificates are carried out in complete detail by hand in
`records/p3/flag-certificates/WORKED-EXAMPLE*.md`.

The changes of variables were found by a systematic search over linear
substitutions and extension fields.  The remaining 480 fields of the census stay
undecided at this date.  The certificates are additive records: later
additions will extend the directory without altering existing files.

The same kind of certificate settles the one field at p = 5 that the
paper left undecided, D = -781922404 (K = Q(sqrt(-195480601)), 5-class
group Z/120 x Z/15 x Z/5): after a change of variables over F_25 the
completed basis of the cubic relation ideal has the four head words
113, 122, 133, 1232, so the 5-class tower group is mild of cohomological
dimension 2, and Theorem D holds for all 204 fields of 5-class rank
three with |D_K| < 2^30.  Over F_5 itself no change of variables
terminates by degree 11 (every one of the 186 flags still acquires a
new head word in degree 11); over F_25 the line of every terminating flag is
the same, namely the non-reduced closed point of degree two of the
norm-degeneracy scheme at which the transversality criterion of the
paper provably fails.  The certificate is
`records/p5/flag-certificates/781922404.json`, its relation tensor is
the entry `cubic_relation_matrix` of `records/p5/D-781922404/result.gp`
(the output of the run whose arithmetic certificate is
`certificates/p5/K-195480601-p5/certificate.gp`), and the same verifier
checks it (`python3 tools/verify_flag_certificate.py
records/p5/flag-certificates/781922404.json`); it was computed by
Singular's Letterplace implementation and recomputed by an independent
completion engine.

The same route also settles the single field of the range with 3-class
rank four, D = -653329427 (class group Z/210 x Z/3 x Z/3 x Z/3), which
the paper reports as undecided because its Bockstein cone is empty: after
a change of variables in GL_4(F_3) the completed basis of the ideal
generated by the four cubic initial forms has twelve head words of degree
at most five, and the words in four letters avoiding them are counted by
1/(1-4z+4z^3), so the 3-class tower group is mild of cohomological
dimension 2 --- to our knowledge the first mild p-class tower group of
p-class rank four for a field in the range of the paper, and the first
one obtained without the transversality criterion (which cannot apply
here, the cone being empty).  Of the 2 080 flags of F_3^4, eight were
found to terminate by degree 9.  The arithmetic behind it is certified:
`records/p3/rank-four-653329427/` now carries the witness certificate
`certificate.gp` built by `tools/build_witness_certificate.gp` and its
check by `tools/verify_certificate_general.gp` (`verification.log`:
812 checks), and `tensor.json` is the tensor that verifier reconstructs,
identical to the operators recorded there in August.  The flag
certificate is `records/p3/flag-certificates-rank-four/653329427.json`;
the rank-generic verifier `tools/verify_flag_certificate_d.py` recomputes
the completion from `tensor.json` (about two minutes in plain Python), and
Singular's Letterplace implementation agrees.

In fact the same combination decides every field of 3-class rank four
known to us.  There are exactly 62 imaginary quadratic fields of
3-class rank four with |D_K| < 2 * 10^10 --- the enumeration by
Belabas's cubic-field counts, which determine the 3-rank
unconditionally through the classical multiplicity (3^r - 1)/2, agrees
field for field, class group for class group, with the unconditional
tables of Mosunov and Jacobson --- and one further field,
D = -20217903567, just beyond that bound.  **Every one of these 63
fields has a mild 3-class tower group of cohomological dimension 2.**
Four carry a transverse point on their Bockstein cone and fall under
the transversality criterion (the three fields announced above on
26 August, and D = -19288632407, which is also decided independently
by a flag); the remaining 59 are settled by flag certificates ---
57 with a completed basis by degree 9, two (D = -6905985272 and
D = -7309564084, and no other) only at degree 11.  In every observed case
the cone is nonempty exactly when the 3-class group has an invariant
divisible by 9 (an element of order 9); for the 35 fields with
elementary 3-class group the cone is always empty.  Each
flag certificate was computed by Singular's Letterplace implementation,
recomputed by an independent plain-Python completion engine, and each
of the 60 flag-decided fields carries the same arithmetic witness
certificate as the field above (812 checks each, base class group
unconditional by `bnfcertify`), with the reconstructed relation tensor
identical to the one the flag search used.  The census, the per-field
records --- relation tensor, flag certificate, witness certificate,
verification log --- and the manifest are under
`records/p3/rank-four-census/`.

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
| `records/p3/flag-certificates/` | one certificate per field decided after the paper: the change of variables, the coefficient field, the head words; format note, worked examples, index |
| `records/p5/flag-certificates/` | the certificate for the field D = -781922404 at p = 5, with its relation tensor under `records/p5/D-781922404/` |
| `records/p3/flag-certificates-rank-four/` | the certificate for the rank-four field D = -653329427; its witness certificate, verification log and tensor are under `records/p3/rank-four-653329427/` |
| `records/p3/rank-four-census/` | the rank-four census through 2 * 10^10 plus the next field D = -20217903567: manifest and, per field, the relation tensor, the flag certificate, the witness certificate and its verification log |
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
