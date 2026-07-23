// MIT License

#ifndef MASSEY_TENSOR_H
#define MASSEY_TENSOR_H

#include <pari/pari.h>

/*
 * Minimal quadratic data:
 *   [m, D_basis, D_pairs]
 * D_basis[i] = D_{e_i}.
 * D_pairs are ordered lexicographically by (i,j), 1 <= i < j <= m,
 * and store D_{e_i+e_j}.
 */
GEN my_reconstruct_secondary_norm(
    GEN p, GEN quadratic_family, GEN t);
GEN my_secondary_norm_delta_basis(
    GEN p, GEN quadratic_family, long i, long k);
GEN my_triple_massey_word_matrix(
    GEN p, GEN quadratic_family);
GEN my_triple_massey_contract(
    GEN p, GEN word_matrix, long m, GEN x, GEN y, GEN z);
void my_validate_triple_massey_identities(
    GEN p, GEN word_matrix, long m);

#endif // MASSEY_TENSOR_H
