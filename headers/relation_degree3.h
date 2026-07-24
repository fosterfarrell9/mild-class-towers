// MIT License

#ifndef RELATION_DEGREE3_H
#define RELATION_DEGREE3_H

#include <pari/pari.h>

/**
 * Analyze the fixed 3-by-27 degree-three relation matrix over F_5.
 *
 * This is pure finite-field/noncommutative algebra: columns are the 27 words
 * in `a,b,c`, while rows are the three relation functionals extracted from
 * the arithmetic tensor.
 */
void my_run_relation_degree3_fixture(GEN T, GEN p);

/**
 * Search and verify the stored Anick strong-freeness certificate.
 *
 * The routine tests admissible word orders and checks the selected leading
 * monomials combinatorially; it performs no number-field arithmetic.
 */
void my_run_mild_certificate_fixture(GEN T, GEN p);

/**
 * Search for and exactly verify an Anick strong-freeness witness.
 *
 * Sparse-first ordered bases of F_p^3 and all six degree-lex variable orders
 * are tested. Ordered row reduction supplies the relation basis change and
 * leaders, which are accepted only by the generic overlap test.
 *
 * `candidate_limit` counts invertible degree-one bases. A negative value
 * requests exhaustive enumeration of GL(3,F_p) within the six implemented
 * degree-lex orders.
 *
 * @return Empty if no configured-search witness is found; otherwise
 * `[M, U, U*T_M, leading_words, variable_order, 1]`.
 */
GEN my_find_strongly_free_witness(GEN T, GEN p, long candidate_limit);

/** Test combinatorial freeness for three encoded length-three words. */
int my_anick_words_combinatorially_free(GEN words);

#endif
