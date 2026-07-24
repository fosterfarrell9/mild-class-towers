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

#ifndef EXT_AND_AUT_H
#define EXT_AND_AUT_H

#include <pari/pari.h>

/**
 * Build relative/absolute models for a list of class-field polynomials.
 *
 * Each output entry is `[Labs, Lrel, sigma, Lbnr]`.  `Labs` and `Lbnr` are
 * relative-field BNF/BNR data used by the search routines; `Lrel` is the
 * compatible relative model; and `sigma` is a nontrivial automorphism fixing
 * the base field.  At this stage sigma is only a generator candidate.  The
 * prescribed-character secondary-norm code later applies the paper's Artin
 * normalization.
 *
 * @param base Base-field BNF/nf accepted by PARI's relative-field routines.
 * @param base_clf Vector of relative class-field defining polynomials.
 * @param p Expected prime relative degree.
 * @param p_rk Number of extension entries to construct.
 * @param D_prime_vect Auxiliary discriminant-prime data for `rnfpolredbest`.
 * @return A stack-independent vector of four-component extension records.
 */
GEN my_ext(GEN base, GEN base_clf, GEN p, int p_rk, GEN D_prime_vect);

#endif // EXT_AND_AUT_H
