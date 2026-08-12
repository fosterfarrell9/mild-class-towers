# Arithmetic certificate for \(K=\mathbf Q(\sqrt{-2800905})\), \(p=5\)

This directory contains the public arithmetic certificate for the explicit
example in *Secondary Norms, Triple Massey Products, and Mild Unramified
\(p\)-Class Tower Groups*. It supports the six secondary-norm matrices for

```text
p = 5
K = Q(s), s^2 + 2800905 = 0
disc(K) = -11203620.
```

This is a machine-checkable certificate based on PARI's exact number-field
arithmetic. It is not a formal proof and does not claim independence from
PARI.

## Files

- `certificate.gp` is a textual, PARI-readable certificate. It contains the
  18 main entries (`a`, `b`, `c`, `a+b`, `a+c`, `b+c`, each at `e1`, `e2`,
  `e3`), like every certificate of the collection.  (Until 2026-08-14 this
  file additionally carried nine doubled-character entries `2a`, `2b`,
  `2c`; they were removed so that all certificates carry exactly the
  entries the paper's certificate definition requires.  The removed
  entries remain in the repository history, and the doubled-character
  audit is documented with the searches in the paper's appendix.)
- `secondary-norms.gp` is a small derived data file with the six verified
  secondary-norm matrices in the format of the `result.gp` records,
  regenerated deterministically by `gp -qf tools/export_secondary_norms.gp`;
  it lets finite-algebra tools read all computed fields uniformly.
- The standalone verifier `verify_certificate.c` and its `Makefile` are
  shared between all certified fields and live one directory up, in
  `certificate/`.

The verifier does not run the search and does not construct a BNF for any
degree-five relative field \(L_x/K\). In particular, successful verification
does not require complete knowledge of \(\mathrm{Cl}(L_x)\) or
\(O_{L_x}^{\times}\).

## Logical separation

The full Massey-pari computation is used only to find and export explicit
objects:

```text
relative-field search
    -> persisted I', compact t_AC, a', J, sigma, and q
    -> standalone exact verifier.
```

The verifier follows this chain:

```text
certified base BNF
    -> stored character vector and released label convention
    -> exact Artin kernel identifies L_x/K
    -> exact Artin exponents identify normalized sigma_x
    -> (1-sigma)^2 I' (t_AC) i(J) = O_L
    -> N(t_AC)/a' is a base-field unit
    -> exact reduction at the stored odd prime proves the sign is +1
    -> N(t_AC) = a'
    -> certified base-field coordinates of [N(I')].
```

The AC1 orientation is the paper's orientation. The same uninverted \(I'\)
used in AC1 is used in the ideal norm.

For each entry the verifier also requires the stored relative and absolute
polynomials to define the same concrete PARI model: the canonical absolute
polynomial inside the reconstructed `rnf` must equal the polynomial used to
construct the absolute `nf`.

The verifier is field-generic.  It accepts certificates with the 18 main
entries (as released throughout the collection) and also tolerates
additional doubled-character entries (27 in total, the historical form
of this file); each present label/column pair is required exactly once,
and the doubled-character checks run only when the `2a`, `2b`, `2c`
entries are present.  The comparison against the published matrices of the
paper is performed only when the certified discriminant is the
principal example's; for any other field the six reconstructed matrices
are printed, and an optional second command-line argument naming a
`result.gp` record cross-checks them against its
`secondary_norm_samples` (`RESULT_RECORD_MATCH=PASS`).  See
`REPRODUCING.md` at the repository root for how to export certificates
for further fields.

## Building and verifying

PARI 2.17.4 is required for the released representation and is checked by the
verifier. As elsewhere in this repository, the PARI build must expose the
exact internal `rnfcycaut`, `allauts`, and `cyclicrelfrob` routines declared in
`headers/pari_internal.h`; the required one-file patch of stock
PARI 2.17.4 and the corresponding build steps are documented in
`doc/pari-2.17.4-patch.md` at the repository root. From a clean checkout:

```sh
cd certificate
make PARI=/path/to/pari-prefix
./verify_certificate K-2800905-p5/certificate.gp
```

The patched PARI can be built with `tools/build-patched-pari.sh`, which
installs it into `$HOME/.local` by default:

```sh
tools/build-patched-pari.sh
cd certificate && make PARI="$HOME/.local"
./verify_certificate K-2800905-p5/certificate.gp
```

Successful output starts with:

```text
BASE_BNF_CERTIFIED=PASS
```

then reports `AC1=PASS`, `AC2=PASS`, and the norm class for all 18 entries,
only after also reporting
`FIELD_MODEL_COMPATIBILITY=PASS`, `ARTIN_CHARACTER=PASS`, and
`SIGMA_NORMALIZATION=PASS`. It is followed by:

```text
D_a=PASS
D_b=PASS
D_c=PASS
D_(a+b)=PASS
D_(a+c)=PASS
D_(b+c)=PASS
D_(2a)=4D_a PASS
D_(2b)=4D_b PASS
D_(2c)=4D_c PASS

CERTIFICATE VERIFIED
```

Any failed assertion names the affected character and column.

## Certificate representation

`certificate.gp` is one GP vector expression, split over lines for
readability. It uses only textual GP encodings; it contains no raw PARI
pointers or process-local memory layouts.

The top-level schema is:

```text
[
  format_version,
  generator_PARI_VERSION_CODE,
  p,
  base_polynomial,
  base_discriminant,
  [class_cyclic_invariants, class_number, class_generators,
   torsion_unit_order, torsion_unit_generator],
  entries
]
```

Each entry is:

```text
[
  character_label,
  column_number,
  prescribed_character_vector,
  relative_polynomial,
  absolute_polynomial,
  normalized_sigma_in_absolute_integral_basis,
  a_prime_in_base_integral_basis,
  J_as_base_ideal_HNF,
  I_prime_as_absolute_ideal_HNF,
  t_AC_as_[factor_column,signed_exponent_column],
  rational_prime_ell,
  prime_ideal_q_above_ell,
  expected_norm_class_column_mod_5
]
```

The base class coordinates are relative to the exact class-group generator
HNFs stored in the header. The verifier recomputes and fully certifies the
base BNF, then requires its cyclic invariants, class number, and generator
HNFs to match that convention.

The released character convention on those three certified \(5\)-relevant
class-group generators is:

```text
a=[1,0,0]    b=[0,1,0]    c=[0,0,1]
a+b=[1,1,0]  a+c=[1,0,1]  b+c=[0,1,1]
2a=[2,0,0]   2b=[0,2,0]   2c=[0,0,2]
```

The verifier computes exact Artin symbols of the three base generators in the
stored extension. Relative to PARI's exact cyclic relative generator, their
exponent vector must span precisely the line of the stored character. The
verifier then identifies the stored sigma as an exact power of that relative
generator and requires the Artin exponent vector relative to stored sigma to
equal the character vector coordinate-for-coordinate. Thus a different
nontrivial power of sigma is rejected even though it still fixes \(K\) and has
order five.

Absolute ideal matrices and basis columns are interpreted using `nfinit` of
the stored absolute polynomial. Relative objects are interpreted using
`rnfinit` of the stored relative polynomial over the certified base field.
The canonical `rnf_get_polabs` polynomial is required to equal the stored
absolute polynomial before either model is used for cross-model conversions.
The compact element is never expanded: its principal ideal is reconstructed
with `famat_idealfactor` and `idealfactorback`, and its residue is evaluated
factor by factor using exact finite-field exponentiation.

## Regenerating the certificate

Regeneration intentionally invokes the full search and is separate from
verification:

```sh
make PARI=/path/to/pari-prefix
MASSEY_ARITHMETIC_AUDIT=1 \
MASSEY_CERTIFICATE_EXPORT=certificate/K-2800905-p5/certificate.gp \
./build/massey 5 's^2+2800905'
```

The exporter writes an entry only after the existing arithmetic audit has
checked AC1, AC2, and the independently extracted norm class. Regeneration is
expensive; verification from the persisted certificate is the intended
public workflow.

The hardening checks use fields already present in format version 1.
`certificate.gp` therefore did not require regeneration.

## Human-readable guide

`certificate-guide.pdf` (source `certificate-guide.tex`) restates the
content of `certificate.gp` entry by entry in the notation of the
paper, including per-entry facts recomputed by the generating script.
It is generated deterministically from the certificate; from the
repository root:

```sh
gp -qf tools/certificate_guide.gp
cd certificate/K-2800905-p5 && pdflatex certificate-guide.tex
```

The generator reads the environment variable `CERT_DIR` to produce the
guide of any other certified field (it defaults to this directory).

The guide is exposition only; the verification chain is implemented
solely by `verify_certificate.c`.
