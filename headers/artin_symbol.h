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

#ifndef ARTIN_SYMBOL_H
#define ARTIN_SYMBOL_H

#include <pari/pari.h>

/**
 * Compute the Artin exponent of a base-field ideal in a cyclic extension.
 *
 * The extension represented by `Lrel` is assumed cyclic, unramified, and of
 * prime degree `p`.  The result identifies its Artin symbol with an exponent
 * in Z/pZ relative to PARI's `rnfcycaut` generator.  The ideal is factored in
 * the certified/search BNF `K`, and exact relative Frobenius exponents are
 * added with their ideal-factorization multiplicities.
 *
 * This routine identifies the Artin character of an already constructed
 * extension; it does not normalize the generator to a prescribed character.
 * The caller performs that separate normalization.
 *
 * @param Labs Absolute BNF model compatible with `Lrel`.
 * @param Lrel Relative cyclic extension over `K`.
 * @param K Base-field BNF.
 * @param I_K Base fractional ideal, normally in HNF.
 * @param p Prime relative degree and modulus for the returned exponent.
 * @return An integer representative of the Artin exponent modulo `p`.
 */
int my_Artin_symbol (GEN Labs, GEN Lrel, GEN K, GEN I_K, int p);

#endif // ARTIN_SYMBOL_H
