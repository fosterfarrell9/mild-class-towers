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

#ifndef FIND_CUP_MATRIX_H
#define FIND_CUP_MATRIX_H

#include <pari/pari.h>

/**
 * Worker wrapper used by PARI's thread/closure interface.
 *
 * The argument bundle contains one extension and the shared base arithmetic;
 * the return value is the relation contribution for worker index `i`.
 */
GEN compute_my_relations(long i, GEN args);

/** Compute relation data in parallel across the selected extensions. */
int my_relations_par(GEN K_ext, GEN K, GEN p, int p_rk, GEN Ja_vect, int r_rk);

/**
 * Compute cup-product relation matrices from Artin evaluations.
 *
 * Matrix columns are indexed by arithmetic inputs `(a_j,J_j)` and the row
 * indexing is derived from ordered pairs of H^1 basis characters.
 */
int my_relations (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk);

/**
 * Compute higher Massey matrices for the selected subgroup/extension data.
 *
 * The long implementation combines ideal lifts, Artin symbols, and the
 * defining-system depth `n`; it prints the research computation and returns a
 * status code rather than a persistent matrix object.
 */
int my_massey_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups, int n);

/** Print the Massey products and relation words used by the main program. */
void my_print_massey(GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups);

#endif // FIND_CUP_MATRIX_H
