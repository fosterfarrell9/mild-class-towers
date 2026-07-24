// MIT License

#ifndef EXAMPLE_PIPELINE_H
#define EXAMPLE_PIPELINE_H

#include <pari/pari.h>

#define MASSEY_EXAMPLE_STATUS_ARITHMETIC_FAILED \
    "ARITHMETIC_COMPUTATION_FAILED"
#define MASSEY_EXAMPLE_STATUS_RANK_LT_3 "RANK_LT_3"
#define MASSEY_EXAMPLE_STATUS_NO_WITNESS \
    "NO_STRONGLY_FREE_BASIS_FOUND"
#define MASSEY_EXAMPLE_STATUS_PROVED \
    "STRONGLY_FREE_BASIS_FOUND"

/**
 * Compute, verify, and persist one compact rank-three example record.
 *
 * The arithmetic phase certifies `K`, constructs the six quadratic samples of
 * the secondary norm, and requires the exact character, sigma, AC1, AC2, and
 * norm-class audit already used by the reference computation. The separate
 * finite-field phase reconstructs the cubic relation matrix and searches for
 * an Anick strong-freeness witness.
 *
 * @return The GP result record copied off the local PARI stack.
 */
GEN my_compute_example_result(
    GEN K, GEN p, GEN polynomial, GEN discriminant,
    GEN Ja_vect, GEN D_prime_vect, const char *output_path,
    long strong_search_limit);

#endif
