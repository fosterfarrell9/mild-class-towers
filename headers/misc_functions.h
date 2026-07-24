// MIT License

// Copyright (c) 2025 [Eric Ahlqvist]

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef MISC_FUNCTIONS_H
#define MISC_FUNCTIONS_H

#include <pari/pari.h>

/** Print PARI's runtime type name for debugging. */
void print_pari_type(GEN x);

/** Print whether the base class group has the requested p-rank. */
void my_test_p_rank (GEN K, int p_int);

/** Abort unless the supplied absolute number field is Galois. */
void my_check_galois(GEN K);

/** Concatenate the rows of two PARI matrices with equal column count. */
GEN concatenate_rows(GEN M1, GEN M2);

/**
 * Apply the ideal operator `1-sigma` to `I`.
 *
 * Returns `I * sigma(I)^(-1)`.  Applying this helper twice gives the
 * `(1-sigma)^2 I'` term in the paper-oriented AC1 equation.
 *
 * @param L Absolute nf/BNF model.
 * @param sigma Integral-basis image defining the automorphism.
 * @param I Fractional ideal, normally represented by an HNF matrix.
 * @return A stack-independent fractional ideal.
 */
GEN my_1MS_ideal (GEN L, GEN sigma, GEN I);

/** Flatten a PARI ideal factorization to `[prime ideals, exponents]`. */
GEN my_find_primes_in_factorization(GEN LyAbs, GEN factorization);

/**
 * Compute the class-group matrix of `(1-sigma)^n`.
 *
 * Rows and columns use the BNR class-group generator convention of `Lbnr`.
 * This relative-BNF computation belongs to candidate search, not to the
 * standalone exact certificate verifier.
 */
GEN my_1MS_operator_2 (GEN Labs, GEN Lbnr, GEN sigma, int n);

/**
 * Apply N_{L/K} factorwise to a compact `nf_GENMAT` element.
 *
 * The signed exponent column is preserved; only each absolute-basis factor is
 * converted to the relative model, normed, and returned in the base integral
 * basis.  The potentially enormous represented element is never expanded.
 */
GEN my_rel_norm_compact(GEN Labs, GEN Lrel, GEN K, GEN compact_elt);

/** Form a linear combination of integer column vectors with exponents `exp`. */
GEN my_vect_from_exp (GEN basis, GEN exp);

/**
 * Solve the class-group equation `(1-sigma)^n[I] = [iJ]`.
 *
 * Returns PARI's modular-solver description `[I_0, kernel]` in exponent
 * coordinates of the relative BNF class generators.  This is a search helper:
 * later code materializes candidates and exact auditing verifies them.
 */
GEN my_H90_2 (GEN L, GEN iJ, GEN oneMS_operator, int n);

/** Return fundamental (and when relevant torsion) unit generators modulo p. */
GEN my_find_units_mod_p (GEN K, GEN p);

/**
 * Matrix of relative norms on chosen extension-unit generators.
 *
 * Columns correspond to fundamental units followed by the torsion generator
 * of `Labs`; rows are base-unit coordinates returned by `bnfisunit`.
 */
GEN my_norm_operator (GEN Labs, GEN Lrel, GEN K, GEN p);

/** Return reduced ideals generating the p-relevant part of Cl(K). */
GEN my_find_p_gens (GEN K, GEN p);

/**
 * Build Ahlqvist--Carlson input pairs `(a',J)`.
 *
 * Each pair satisfies `div(a') + pJ = 0`, equivalently `(a')J^p=O_K`.
 * Unit inputs are represented with `J=O_K`.
 */
GEN my_find_Ja_vect(GEN K, GEN J_vect, GEN p, GEN units_mod_p);

/** Enumerate exponent vectors for a finite abelian group with invariants `cyc`. */
GEN my_get_vect (int n, GEN cyc);

/** Enumerate F_p-linear combinations of the supplied kernel basis. */
GEN my_get_sums (GEN basis, int p);

/**
 * Search for the ideal `I'` in the Ahlqvist--Carlson equations.
 *
 * For each `(a',J)`, relative class and unit data are searched for an ideal
 * whose `(1-sigma)^n` image has the required class and whose compact
 * principal generator can be corrected to the prescribed norm.  For `n=2`,
 * the returned uninverted ideal is the `I'` in
 * `(1-sigma)^2 I' + div(t_AC) + i(J)=0`.
 *
 * This routine finds candidates using `Cl(L)` and `O_L^x`; exact AC1/AC2
 * auditing and the standalone certificate verifier are the proof-critical
 * validation paths.
 */
GEN my_H90_vect_2 (GEN Labs, GEN Lrel, GEN Lbnr, GEN K, GEN sigma, GEN Ja_vect, GEN p, int n);

/** Enumerate all ideal representatives of Cl(K); retained for legacy tests. */
GEN my_get_clgp (GEN K);

/** Print the unramified p-extensions selected from the base class group. */
void my_unramified_p_extensions(GEN K, GEN p, GEN D_prime_vect);

/** Compute base-class lifts associated with an extension's transfer kernel. */
GEN my_ideal_lifts (GEN Labs, GEN Lrel, GEN K, GEN p);

/** Explore unramified p-extensions together with class-group transfer data. */
void my_unramified_p_extensions_with_transfer(GEN K, GEN p, GEN D_prime_vect);

/**
 * Select index-p subgroups whose class fields have the smallest class numbers.
 *
 * This heuristic search constructs the candidate class fields and returns
 * `p_rank` subgroup matrices; it is not a proof-critical certification step.
 */
GEN my_best_subgroups(GEN K, long p_rank, GEN subgroups, GEN D_prime_vect);

#endif // MISC_FUNCTIONS_H
