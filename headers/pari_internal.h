#ifndef MASSEY_PARI_INTERNAL_H
#define MASSEY_PARI_INTERNAL_H

#include <pari/pari.h>

/*
 * Internal PARI functions used by the original Ahlqvist code.
 *
 * These are static in stock PARI 2.17.4 and therefore require
 * a locally patched libpari in which the corresponding `static`
 * qualifiers have been removed.
 */
GEN rnfcycaut(GEN rnf);
GEN allauts(GEN rnf, GEN aut);
long cyclicrelfrob(GEN rnf, GEN auts, GEN pr);

#endif