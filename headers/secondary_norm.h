// MIT License

#ifndef SECONDARY_NORM_H
#define SECONDARY_NORM_H

#include <pari/pari.h>

/* Count cyclic invariant factors divisible by p, directly or from bnf K. */
long my_p_rank_from_cyc(GEN cyc, GEN p);
long my_p_class_rank(GEN K, GEN p);
GEN my_p_relevant_coordinates(GEN K, GEN values, GEN p);
GEN my_subgroup_character(GEN K, GEN H, GEN p);
long my_sigma_exponent(GEN Labs, GEN Lrel, GEN sigma_H90, GEN p);
GEN my_automorphism_power_checked(
    GEN Labs, GEN Lrel, GEN K, GEN sigma, long exponent, GEN p);

/*
 * Return the matrix of D_t : E -> Cl(K)/p.
 * Rows are the p-relevant class-group coordinates and column j is D_t(e_j).
 * This implementation currently supports p > 3 only.
 */
GEN my_secondary_norm_operator(
    GEN K, GEN p, GEN prescribed_character_t,
    GEN Ja_vect, GEN D_prime_vect);

/*
 * Return [m, D_basis, D_pairs], where D_pairs are ordered lexicographically
 * by (i,j), 1 <= i < j <= m.
 */
GEN my_secondary_norm_basis_family(
    GEN K, GEN p, GEN Ja_vect, GEN D_prime_vect);

#endif // SECONDARY_NORM_H
