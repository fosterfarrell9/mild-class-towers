// MIT License

#ifndef MASSEY_TENSOR_H
#define MASSEY_TENSOR_H

#include <pari/pari.h>

/**
 * Reconstruct D_t from its quadratic basis and pair values over F_p.
 *
 * `quadratic_family` is `[m, D_basis, D_pairs]`, where `D_basis[i]` is
 * `D_{e_i}` and `D_pairs` contains `D_{e_i+e_j}` in lexicographic `(i,j)`
 * order for `1 <= i < j <= m`.
 *
 * The returned matrix keeps the convention that columns are inputs `e_j` and
 * rows are p-relevant class coordinates.
 */
GEN my_reconstruct_secondary_norm(
    GEN p, GEN quadratic_family, GEN t);

/** Return the polarized matrix Delta D(e_i,e_k) used in triple products. */
GEN my_secondary_norm_delta_basis(
    GEN p, GEN quadratic_family, long i, long k);

/**
 * Assemble coefficients of degree-three words from the polarized norm.
 *
 * Column `(i,j,k)` in lexicographic word order stores the evaluations of the
 * corresponding triple Massey coefficient on all arithmetic inputs.
 */
GEN my_triple_massey_word_matrix(
    GEN p, GEN quadratic_family);

/** Contract a degree-three word matrix with characters x, y, and z. */
GEN my_triple_massey_contract(
    GEN p, GEN word_matrix, long m, GEN x, GEN y, GEN z);

/** Check outer symmetry, cyclic shuffle, and diagonal tensor identities. */
void my_validate_triple_massey_identities(
    GEN p, GEN word_matrix, long m);

#endif // MASSEY_TENSOR_H
