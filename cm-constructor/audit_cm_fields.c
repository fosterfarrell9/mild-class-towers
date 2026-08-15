#define _POSIX_C_SOURCE 200809L

/* Exact, field-only audit for the experiment's CM output.
 *
 * The Artin-kernel calculation is intentionally the same one used by
 * verifier/verify_certificate.c: rnfcycaut + allauts + cyclicrelfrob on
 * the fixed base-class generators.  No relative class or unit group is used.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "headers/pari_internal.h"

enum {
    CM_FORMAT = 1,
    CM_PARI_VERSION,
    CM_P,
    CM_BASE_POLYNOMIAL,
    CM_DISCRIMINANT,
    CM_BASE_DATA,
    CM_INVARIANT_DATA,
    CM_PRECISION_DATA,
    CM_ENTRIES,
    CM_EVALUATION_SECONDS,
    CM_TOTAL_SECONDS
};

enum {
    CM_ENTRY_LABEL = 1,
    CM_ENTRY_CHARACTER,
    CM_ENTRY_RAW_POLYNOMIAL,
    CM_ENTRY_REDUCED_POLYNOMIAL,
    CM_ENTRY_ROUND_EXPONENT,
    CM_ENTRY_MAX_REAL_DISTANCE,
    CM_ENTRY_MAX_IMAGINARY_DISTANCE,
    CM_ENTRY_COSET_COUNTS,
    CM_ENTRY_EXTRACTION_SECONDS
};

static double
monotonic_seconds(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec / 1.0e9;
}

static void
fail(const char *label, const char *message)
{
    if (label)
        pari_fprintf(stderr, "CM AUDIT FAILURE %s: %s\n", label, message);
    else
        pari_fprintf(stderr, "CM AUDIT FAILURE: %s\n", message);
    pari_close();
    exit(EXIT_FAILURE);
}

static GEN
read_expression(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (!file) fail(NULL, "cannot open input");
    if (fseek(file, 0, SEEK_END) != 0)
        fail(NULL, "cannot seek input");
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0)
        fail(NULL, "cannot size input");
    char *text = malloc((size_t)size + 1);
    if (!text) fail(NULL, "cannot allocate input buffer");
    if (fread(text, 1, (size_t)size, file) != (size_t)size)
        fail(NULL, "cannot read complete input");
    fclose(file);
    text[size] = '\0';
    GEN value = gp_read_str(text);
    free(text);
    return value;
}

static GEN
canonical_artin_exponents(GEN Lrel, GEN K, GEN p, const char *label)
{
    /* Match the verifier's operation order: rnfeltup0 initializes the
     * absolute-field cache in the rnf object before rnfcycaut reads it. */
    (void)rnfeltup0(
        Lrel, pol_x(nf_get_varn(bnf_get_nf(K))), 1);
    GEN relative_automorphism = rnfcycaut(Lrel);
    if (!relative_automorphism)
        fail(label, "rnfcycaut did not return a generator");
    GEN all_automorphisms = allauts(Lrel, relative_automorphism);
    GEN generators = bnf_get_gen(K);
    GEN cyc = bnf_get_cyc(K);
    GEN result = cgetg(4, t_COL);
    long coordinate = 1;

    for (long i = 1; i < lg(cyc); ++i)
    {
        if (!dvdii(gel(cyc, i), p)) continue;
        GEN factorization = idealfactor(K, gel(generators, i));
        GEN exponent = gen_0;
        for (long factor = 1; factor < lg(gel(factorization, 1)); ++factor)
        {
            long frobenius = cyclicrelfrob(
                Lrel, all_automorphisms,
                gmael(factorization, 1, factor));
            exponent = Fp_add(
                exponent,
                Fp_mul(
                    stoi(frobenius),
                    gmael(factorization, 2, factor), p),
                p);
        }
        gel(result, coordinate++) = exponent;
    }
    if (coordinate != 4)
        fail(label, "base p-relevant class rank is not three");
    return result;
}

static void
verify_character_line(GEN canonical, GEN character, GEN p, const char *label)
{
    long pivot = 0;
    for (long i = 1; i <= 3; ++i)
        if (signe(gel(character, i)))
        {
            pivot = i;
            break;
        }
    if (!pivot) fail(label, "zero character");
    GEN scale = Fp_div(gel(canonical, pivot), gel(character, pivot), p);
    if (!signe(scale)) fail(label, "zero Artin character");
    for (long i = 1; i <= 3; ++i)
        if (!equalii(
                gel(canonical, i),
                Fp_mul(scale, gel(character, i), p)))
            fail(label, "exact Artin character does not have ker(x)");
}

static GEN
stored_absolute_polynomial(GEN certificate, const char *label)
{
    GEN entries = gel(certificate, 7);
    for (long i = 1; i < lg(entries); ++i)
        if (strcmp(GSTR(gmael(entries, i, 1)), label) == 0)
            return gmael(entries, i, 5);
    fail(label, "label not found in stored certificate");
    return NULL;
}

int
main(int argc, char **argv)
{
    if (argc < 2 || argc > 3)
    {
        fprintf(
            stderr,
            "usage: %s <cm-output.gp> [stored-certificate.gp]\n",
            argv[0]);
        return 2;
    }

    pari_init_opts(
        1L << 30, 1048576,
        INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    paristack_setsize(1L << 30, 1L << 33);
    double total_started = monotonic_seconds();

    GEN cm = read_expression(argv[1]);
    GEN certificate = argc == 3 ? read_expression(argv[2]) : NULL;
    if (typ(cm) != t_VEC || lg(cm) != 12 || !equaliu(gel(cm, CM_FORMAT), 1))
        fail(NULL, "invalid CM output schema");
    if (!equaliu(gel(cm, CM_PARI_VERSION), PARI_VERSION_CODE))
        fail(NULL, "CM output has a different PARI version");
    if (certificate)
    {
        if (typ(certificate) != t_VEC || lg(certificate) != 8)
            fail(NULL, "invalid stored certificate schema");
        if (!gequal(gel(cm, CM_BASE_POLYNOMIAL), gel(certificate, 4))
            || !gequal(gel(cm, CM_DISCRIMINANT), gel(certificate, 5)))
            fail(NULL, "CM/certificate base fields differ");
    }

    GEN p = gel(cm, CM_P);
    GEN K = Buchall(gel(cm, CM_BASE_POLYNOMIAL), nf_FORCE, DEFAULTPREC);
    if (bnfcertify0(K, 0) != 1)
        fail(NULL, "bnfcertify(K) failed");
    if (!gequal(bnf_get_cyc(K), gmael(cm, CM_BASE_DATA, 1))
        || !gequal(bnf_get_no(K), gmael(cm, CM_BASE_DATA, 2))
        || !gequal(bnf_get_gen(K), gmael(cm, CM_BASE_DATA, 3)))
        fail(NULL, "CM base class-group convention changed");
    if (certificate
        && (!gequal(bnf_get_cyc(K), gmael(certificate, 6, 1))
            || !gequal(bnf_get_no(K), gmael(certificate, 6, 2))
            || !gequal(bnf_get_gen(K), gmael(certificate, 6, 3))))
        fail(NULL, "stored base class-group convention changed");

    pari_printf("BASE_BNF_CERTIFIED=PASS\n");
    GEN entries = gel(cm, CM_ENTRIES);
    if (typ(entries) != t_VEC || glength(entries) < 1
        || glength(entries) > 6)
        fail(NULL, "CM output must contain one to six characters");
    long verified = 0;

    for (long i = 1; i < lg(entries); ++i)
    {
        double started = monotonic_seconds();
        GEN entry = gel(entries, i);
        const char *label = GSTR(gel(entry, CM_ENTRY_LABEL));
        GEN character = gel(entry, CM_ENTRY_CHARACTER);
        GEN polynomial = gel(entry, CM_ENTRY_REDUCED_POLYNOMIAL);
        /* This exact relative reduction is cosmetic and occurs only after
         * the floating-point trace polynomial has been rounded to Z[y]. */
        GEN relative_polynomial = liftall(
            rnfpolredbest(bnf_get_nf(K), polynomial, 0));
        GEN Lrel = rnfinit(bnf_get_nf(K), relative_polynomial);
        GEN absolute = rnf_get_polabs(Lrel);
        GEN Labs = nfinit0(absolute, 0, DEFAULTPREC);

        if (nf_get_degree(Labs) != 10)
            fail(label, "CM compositum does not have absolute degree ten");
        if (!equalii(
                nf_get_disc(Labs),
                powiu(gel(cm, CM_DISCRIMINANT), 5)))
            fail(label, "CM compositum is not everywhere unramified over K");

        GEN canonical = canonical_artin_exponents(Lrel, K, p, label);
        verify_character_line(canonical, character, p, label);

        if (certificate)
        {
            GEN stored = stored_absolute_polynomial(certificate, label);
            GEN isomorphisms = nfisisom(absolute, stored);
            if (gequal0(isomorphisms))
                fail(label, "CM field is not the stored exact degree-ten field");
            pari_printf(
                "%s  DEGREE=PASS DISC=D^5=PASS ARTIN_CHARACTER=PASS "
                "CANONICAL_ARTIN=%Ps EXACT_SAME_EXTENSION_NFISISOM=PASS "
                "EMBEDDINGS=%ld RELATIVE_POLYNOMIAL=%Ps SECONDS=%.9f\n",
                label, gtovec(canonical), glength(isomorphisms),
                relative_polynomial, monotonic_seconds() - started);
        }
        else
            pari_printf(
                "%s  DEGREE=PASS DISC=D^5=PASS ARTIN_CHARACTER=PASS "
                "CANONICAL_ARTIN=%Ps RELATIVE_POLYNOMIAL=%Ps SECONDS=%.9f\n",
                label, gtovec(canonical), relative_polynomial,
                monotonic_seconds() - started);
        ++verified;
    }
    pari_printf(
        "CM FIELD AUDIT %ld/%ld VERIFIED TOTAL_SECONDS=%.9f\n",
        verified, glength(entries), monotonic_seconds() - total_started);
    pari_close();
    return 0;
}
