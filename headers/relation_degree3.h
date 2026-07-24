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

#endif
