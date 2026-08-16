# The arithmetic certificate format for p = 3

## Scope

This document specifies the standalone arithmetic certificate used for
imaginary quadratic fields `K` whose 3-class rank is three.  A certificate is
one textual GP vector expression.  The verifier reads that expression as data;
it does not execute the file and it does not construct a BNF for an auxiliary
field.

The format inherits the integral-basis and exact-audit discipline of the p = 5
format.  Its p = 3 specialization differs in four mathematical places: the
divisor `J` contributes to the accepted class, the relative extensions are
cubic, all restricted-cube coefficients are retained, and six character
evaluations plus the shuffle identities are mandatory.

## Conventions and orientation

Let

```text
C = Cl(K)/3,
E = Cl(K)[3],
V = C^vee,
```

with ordered bases `chi_1, chi_2, chi_3` of `V` and `e_1, e_2, e_3`
of `E`.  The class-group coordinates used by the certificate are PARI's three
3-relevant class-group generators, in their stored order.  They are unrelated
to the cohomological comparison maps customarily denoted by `kappa_i`.

The unadorned Massey product is the Ahlqvist--Carlson convention
`M = M^AC`.  The global convention comparison is

```text
M^AC(x,y,z) = -M^V(x,y,z).
```

The certificate fixes the relation-row orientation by declaring

```text
bar(r)_ell = - trg^vee(lambda_ell),
lambda_ell(alpha) = <kappa_2(alpha), e_ell>_AV.
```

Consequently its stored expected tensor has

```text
T[ell,i,j,k] = epsilon_{ijk,3}(r_ell)
             = <M^AC(chi_i,chi_j,chi_k), e_ell>_AV.
```

This declaration fixes the sign of the relation classes; it is not an
instruction to invert already oriented relators.

For a nonzero character `x`, the generator `sigma_x` of the unramified cyclic
cubic extension `L_x/K` is oriented by

```text
Art_{L_x/K}(c) = sigma_x ^ x(c mod 3).
```

The verifier recomputes this equality on the three stored base class-group
generators.  A different nontrivial power of `sigma_x` is rejected.

## The exact p = 3 entry identities

An input `e in E` is represented by `(a_prime,J)`, where `a_prime` is a
base-field element and `J` is a fractional ideal/divisor.  An entry stores
`t in L_x^*` and a fractional ideal/divisor `I_prime` of `L_x`.  Additively,
the required identities are

```text
div(a_prime) + 3 J = 0,
(1-sigma_x)^2 I_prime + div(t) + i_x(J) = 0,       (AC1)
N_{L_x/K}(t) = a_prime.                            (AC2)
```

Equivalently, in multiplicative ideal notation,

```text
(a_prime) J^3 = O_K,
(1-sigma_x)^2 I_prime (t) i_x(J) = O_L,
N_{L_x/K}(t) = a_prime.
```

Here `(1-sigma)I` means `I / sigma(I)`, applied twice.  The stored compact
factorization represents the Ahlqvist--Carlson element `t` itself.  The search
oracle initially produces the inverse convention and the builder normalizes
it before emission.

The p = 3 Ahlqvist--Carlson evaluation is

```text
<M(x,x,y),(a_prime,J)>_AV = <y, N_{L_x/K}(I_prime) + J>.
```

Therefore the certified secondary norm value is

```text
D_x(e) = [N_{L_x/K}(I_prime) + J]
       = [N_{L_x/K}(I_prime) J] in Cl(K)/3.
```

The class of `N(I_prime)` alone is not an accepted p = 3 value.  AC2 is
checked twice: its quotient is required to generate the unit ideal, and exact
reduction at a stored odd prime ideal must give `+1`, not `-1`.  Since these
rank-three imaginary quadratic base fields have unit group `{+1,-1}`, this
fixes the unit exactly.

## Top-level schema

The GP vector is

```text
[
  2,
  generator_PARI_VERSION_CODE,
  3,
  base_polynomial,
  base_discriminant,
  [
    class_cyclic_invariants,
    class_number,
    class_generators,
    torsion_unit_order,
    torsion_unit_generator,
    base_integral_basis,
    expected_tensor_3_by_27
  ],
  entries
]
```

The expected tensor rows use `e_1,e_2,e_3`.  Its 27 columns use the ordered
words `X_i X_j X_k` in lexicographic order with `k` fastest.  The orientation
is the negative transgression-dual convention above.

The certificate contains exactly 18 entries: three input columns for each of

```text
x1, x2, x3, x1+x2+x3, x1+x2, x1+x3.
```

No label/column pair may be missing or duplicated.  Certificates whose
format identifier differs from the expected value, additional
character families, and underspecified vectors are rejected.

## Entry schema

Each entry is

```text
[
  character_label,
  input_column,
  prescribed_character_vector,
  relative_polynomial,
  absolute_polynomial,
  normalized_sigma_in_absolute_integral_basis,
  a_prime_in_base_integral_basis,
  J_as_base_ideal_HNF,
  I_prime_as_absolute_ideal_HNF,
  t_as_[factor_column,signed_exponent_column],
  rational_prime_ell,
  prime_ideal_q_above_ell,
  expected_J_corrected_class_mod_3,
  absolute_integral_basis
]
```

Elements of an explicitly stored integral basis are rational numbers or
polynomials in the field generator, of degree smaller than the field degree.
For every stored basis, the verifier expresses its elements in the basis
chosen by its own `nfinit`, requires an integral change-of-basis matrix, and
requires determinant `+1` or `-1`.  Extension-field coordinate vectors,
ideals, and compact factors are converted through that matrix before use.
The quadratic base-field basis receives the same integral and unimodular
audit and is additionally required to agree with the local basis, as in the
p = 5 verifier, because the stored PARI prime-ideal records and base
class-group convention are basis-dependent.  No LLL-reduced basis is assumed
to be canonical.

## Six-character reconstruction and shuffle lock

Write `D1,D2,D3,D123,D12,D13` for the six verified 3 by 3 matrices; rows are
class coordinates and columns are the inputs `e_ell`.  Set

```text
B12 = D1 + D2 - D12,
B13 = D1 + D3 - D13,
B23 = D1 + D2 + D3 - B12 - B13 - D123.
```

For every relation row `ell` and middle index `j`, reconstruction uses

```text
T[ell,i,j,i] = Di[j,ell],
T[ell,i,j,k] = T[ell,k,j,i] = Bik[j,ell]  for i < k.
```

The verifier then checks all ordered triples against both shuffle identities

```text
M(x,y,z) = M(z,y,x),
M(x,y,z) + M(y,z,x) + M(z,x,y) = 0.
```

In characteristic three these imply `M(x,y,x)=M(x,x,y)`.  Equivalently,
for the polarization `Delta D(x,z)=D_x+D_z-D_{x+z}`,

```text
Delta D(x,x) = 2D_x - D_{2x} = D_x,
```

not `-2D_x` as for primes greater than three.  The p = 5 doubled-character
redundancy degenerates because `2x=-x`; the six-character shuffle audit is
therefore mandatory.  Finally, the reconstructed 81 entries must equal the
stored expected tensor entry by entry.

## Acceptance boundary

The relative `bnfinit(L_x)` used by the builder is an oracle only.  The
standalone verifier accepts an entry only after checking:

1. the format identifier, prime 3, PARI representation version, and complete schemas;
2. `bnfcertify(K)`, the base discriminant, class group, unit data, and basis;
3. equality of the relative and absolute defining polynomials, relative degree three,
   and `disc(L_x)=disc(K)^3`;
4. exact Artin kernel identification and exact `sigma_x` normalization;
5. the two integral-basis audits;
6. `(a_prime)J^3=O_K`, AC1, and AC2 including its sign;
7. independent recomputation of `[N(I_prime)J]` in the certified base class
   group;
8. completeness of the 18 values, reconstruction, every shuffle identity,
   and entrywise equality with the expected tensor.

No auxiliary class-group value from the search is read or trusted by the
verifier.  A failed search can leave a field undecided; it cannot produce a
negative result.
