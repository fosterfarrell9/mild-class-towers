// MIT License

#ifndef SECONDARY_NORM_H
#define SECONDARY_NORM_H

#include <pari/pari.h>

/** Count cyclic invariant factors contributing to the p-relevant quotient. */
long my_p_rank_from_cyc(GEN cyc, GEN p);

/** Return the number of cyclic factors of Cl(K) divisible by `p`. */
long my_p_class_rank(GEN K, GEN p);

/**
 * Project class exponent data to the p-relevant coordinates modulo `p`.
 *
 * Coordinates retain the order of `bnf_get_cyc(K)` after omitting invariant
 * factors not divisible by `p`.  The result is a column vector.
 */
GEN my_p_relevant_coordinates(GEN K, GEN values, GEN p);

/**
 * Recover the normalized F_p character whose kernel is the subgroup `H`.
 *
 * The coordinates use the p-relevant certified/search class-group generators
 * of `K`.  The one-dimensional nullspace is scaled so its first nonzero entry
 * is one.
 */
GEN my_subgroup_character(GEN K, GEN H, GEN p);

/**
 * Express an absolute automorphism as a power of PARI's relative Artin generator.
 *
 * @return The unique exponent in `1..p-1`; raises a PARI error if no such
 * power exists.
 */
long my_sigma_exponent(GEN Labs, GEN Lrel, GEN sigma_H90, GEN p);

/**
 * Form and structurally validate a nonzero power of a cyclic automorphism.
 *
 * The returned integral-basis image fixes `K`, has order `p`, and generates
 * Gal(L/K).  The result is copied off the function's local PARI stack.
 */
GEN my_automorphism_power_checked(
    GEN Labs, GEN Lrel, GEN K, GEN sigma, long exponent, GEN p);

/**
 * Compute the prescribed-character secondary norm D_t.
 *
 * The search constructs the unramified cyclic extension cut out by the line
 * F_p t, normalizes sigma so `Art(c) = sigma^{t(c)}`, and searches relative
 * BNF class/unit data for Ahlqvist--Carlson candidates `(I',t_AC)`.  When the
 * arithmetic audit is enabled, those candidates are separately checked by
 * exact AC1/AC2 arithmetic before their norm classes are accepted.
 *
 * Return the matrix of `D_t : E -> Cl(K)/p`.
 * Rows are the p-relevant class-group coordinates and column j is D_t(e_j).
 * This implementation currently supports p > 3 only.
 *
 * @param prescribed_character_t Coordinates on the p-relevant class-group
 * generators, reduced modulo `p` internally.
 * @param Ja_vect Ahlqvist--Carlson pairs `(a',J)` with div(a') + pJ = 0.
 * @return A PARI matrix copied off the local stack.
 */
GEN my_secondary_norm_operator(
    GEN K, GEN p, GEN prescribed_character_t,
    GEN Ja_vect, GEN D_prime_vect);

/**
 * Compute the minimal quadratic family determining all secondary norms.
 *
 * Return [m, D_basis, D_pairs], where D_pairs are ordered lexicographically
 * by `(i,j)`, `1 <= i < j <= m`.  `D_basis[i]=D_{e_i}` and each pair entry is
 * `D_{e_i+e_j}`.  Characters scalar-related over F_p define the same extension
 * kernel but use the correspondingly normalized automorphism.
 */
GEN my_secondary_norm_basis_family(
    GEN K, GEN p, GEN Ja_vect, GEN D_prime_vect);

/**
 * Require exact auditing for subsequently computed secondary norms.
 *
 * The compact example pipeline uses this process-local switch so it cannot
 * accept values produced only by relative-BNF candidate search.
 */
void my_secondary_norm_require_exact_audit(int required);

#endif // SECONDARY_NORM_H
