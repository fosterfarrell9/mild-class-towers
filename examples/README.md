# Compact example records

This directory is reserved for compact, machine-readable computations of
additional imaginary quadratic examples. The detailed release in
`certificate/K-2800905-p5/` remains the reference certificate: compact records
do not duplicate its large `I'` and `t_AC` witnesses.

Compact does not mean weaker arithmetic. Before any secondary-norm matrix is
recorded, the program certifies the base BNF and checks the character kernel,
the prescribed Artin normalization of sigma, AC1 and AC2 by exact ideal
arithmetic, and the independently recovered norm class. Relative BNF/BNR data
are used only to find candidates.

Run one field from the repository root:

```sh
mkdir -p examples/p5/K-2800905
./build/massey --example-result \
  examples/p5/K-2800905/result.gp 5 's^2+2800905'
```

The program does not create parent directories. A result is written only after
the arithmetic and finite-field stages complete.

## Audited example at discriminant -18397407

The directory `p5/K-18397407/` records the computation for

```text
p = 5
K = Q[s]/(s^2-s+4599352)
disc(K) = -18397407
Cl(K) = [40,10,5].
```

It was produced with:

```sh
./build/massey --example-result \
  examples/p5/K-18397407/result.gp 5 's^2-s+4599352'
```

Every one of the 18 secondary-norm columns passed the prescribed-character,
normalized-Artin-generator, exact AC1/AC2, and independent norm-class audit.
The resulting cubic relation matrix has rank 3.  The generic search found a
verified strongly-free witness with leading words `["bba","bcc","bca"]`,
proving `MILD="PROVED"` and `CD=2`.

`result.gp` is the compact machine-readable record.  For this example,
`run.log` is also retained as the complete console transcript of the exact
arithmetic audit.

The strong-freeness search tests 250,000 invertible degree-one bases by
default. Set another positive bound with:

```sh
./build/massey --example-result result.gp \
  --strong-search-limit 500000 5 's^2+2800905'
```

Request every element of `GL_3(F_p)` with:

```sh
./build/massey --example-result result.gp \
  --strong-search-limit exhaustive 5 's^2+2800905'
```

For `p=5`, `|GL_3(F_5)| = 1,488,000`, so the latter tests all 1,488,000
degree-one basis changes against the six implemented degree-lex variable
orders. Even this exhaustive GL search does not prove that no strongly-free
presentation exists outside the implemented witness class.

## GP schema

`result.gp` is a GP vector of `[name,value]` pairs with format version 1:

- `status`: `ARITHMETIC_COMPUTATION_FAILED`, `RANK_LT_3`,
  `NO_STRONGLY_FREE_BASIS_FOUND`, or `STRONGLY_FREE_BASIS_FOUND`;
- `p`, base polynomial and discriminant;
- certified class-group invariants, class number, and PARI generators;
- the character basis as columns;
- the six matrices in order
  `D_a,D_b,D_c,D_(a+b),D_(a+c),D_(b+c)`;
- reconstructed `D_(2a),D_(2b),D_(2c)` scaling checks;
- `arithmetic_exact_audit`;
- the 3-by-27 cubic relation matrix, whose columns are `X_i X_j X_k` with
  the third index varying fastest;
- its rank;
- either an empty strong-freeness witness or
  `[M,U,U*T_M,leading_words,variable_order,1]`;
- the candidate limit and whether it exhausts `GL_3(F_p)`;
- `MILD`, equal to `"PROVED"` only when the recorded Anick witness verifies
  and `"UNKNOWN"` otherwise;
- `CD`, equal to 2 for a verified mildness witness and `"UNKNOWN"` otherwise.

Here the columns of `M` give the images of the old degree-one generators in
the new basis, and `U` changes the relation basis. The search enumerates
sparse ordered bases first and tests all six degree-lex variable orders.

`MILD="PROVED"` means the cubic initial forms have combinatorially free leading
monomials, so Anick's criterion proves strong freeness and the tower group is
mild with cohomological dimension 2. `MILD="UNKNOWN"` and `CD="UNKNOWN"` mean
only that this implemented sufficient-criterion search found no witness; they
are not negative mathematical conclusions.

The current single-field command exits nonzero and need not write a record when
arithmetic fails. A future batch driver will catch that failure and persist
`status="ARITHMETIC_COMPUTATION_FAILED"`; the arithmetic routine itself does
not currently serialize failure records.
