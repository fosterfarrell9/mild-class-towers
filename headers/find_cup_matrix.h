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

//-----------------------------------------------------------------------------

//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>

//------------------------------------------------
// Parallel version of my_relations
//------------------------------------------------

/*
GP;install("compute_my_relations","Gp","./stdin.so");
*/

// Wrapper function for parallel computation
GEN compute_my_relations(long i, GEN args);

int my_relations_par(GEN K_ext, GEN K, GEN p, int p_rk, GEN Ja_vect, int r_rk);

int my_relations (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk);

int my_massey_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups, int n);

void my_print_massey(GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups);

#endif // FIND_CUP_MATRIX_H
