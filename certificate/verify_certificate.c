/**
 * @file verify_certificate.c
 * @brief Standalone exact verifier for the published arithmetic certificate.
 *
 * The verifier reconstructs each character kernel and stored cyclic extension,
 * normalizes sigma to the prescribed character via exact Artin data, reconciles
 * the relative and absolute models, checks AC1 and AC2, and independently
 * recovers the asserted norm class.  It deliberately performs no relative
 * BNF/BNR candidate search: every accepted entry follows the complete exact
 * verification chain encoded in the certificate.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../headers/pari_internal.h"

enum {
    CERT_FORMAT = 1,
    CERT_PARI_VERSION,
    CERT_P,
    CERT_BASE_POLYNOMIAL,
    CERT_DISCRIMINANT,
    CERT_BASE_DATA,
    CERT_ENTRIES
};

enum {
    ENTRY_CHARACTER = 1,
    ENTRY_COLUMN,
    ENTRY_CHARACTER_VECTOR,
    ENTRY_RELATIVE_POLYNOMIAL,
    ENTRY_ABSOLUTE_POLYNOMIAL,
    ENTRY_SIGMA,
    ENTRY_A_PRIME,
    ENTRY_J,
    ENTRY_I_PRIME,
    ENTRY_T_AC,
    ENTRY_ELL,
    ENTRY_PRIME,
    ENTRY_NORM_CLASS,
    ENTRY_INTEGRAL_BASIS       /* format 2 only */
};

static void
fail(const char *label, long column, const char *assertion)
{
    if (label)
        pari_fprintf(
            stderr, "CERTIFICATE FAILURE %s/e%ld: %s\n",
            label, column, assertion);
    else
        pari_fprintf(stderr, "CERTIFICATE FAILURE: %s\n", assertion);
    pari_close();
    exit(EXIT_FAILURE);
}

/*
 * Certificates record elements and ideals as coordinates with respect to an
 * integral basis -- and PARI's integral basis is LLL-reduced, hence not
 * canonical.  nfinit on another machine may hand back a different basis, and
 * the same coordinate vector then denotes a different algebraic number: the
 * stored automorphism stops fixing K, the stored ideals become other ideals,
 * and the verification fails on data that is perfectly correct.  Measured on
 * two machines, eleven of twenty-one certificates were affected.
 *
 * Format 2 therefore records the basis that was used, and the verifier
 * translates the stored coordinates into its own.  The columns of the matrix
 * below are the stored basis elements expressed in the local basis.  Insisting
 * that it be integral with determinant +-1 is what makes this safe: the two
 * bases then span the same ring of integers, so a certificate cannot smuggle
 * in a different order by declaring a convenient "basis".
 */
static GEN
basis_change_matrix(
    GEN nf, GEN stored_basis, const char *label, long column)
{
    long degree = nf_get_degree(nf);
    if (typ(stored_basis) != t_VEC || lg(stored_basis) != degree + 1)
        fail(label, column, "stored integral basis has the wrong length");

    GEN matrix = cgetg(degree + 1, t_MAT);
    for (long i = 1; i <= degree; ++i)
        gel(matrix, i) = algtobasis(nf, gel(stored_basis, i));
    if (!RgM_is_ZM(matrix))
        fail(label, column,
             "stored integral basis is not integral in the local basis");
    if (!is_pm1(ZM_det(matrix)))
        fail(label, column, "stored integral basis is not unimodular");
    return matrix;
}

/**
 * An element given by stored coordinates, in local coordinates.  Coordinates
 * of algebraic numbers are rational in general -- a' carries denominators of
 * the order of 2^40 -- so this must not use the integer-only routines.
 */
static GEN
transformed_element(GEN matrix, GEN coordinates)
{
    return typ(coordinates) == t_COL
        ? RgM_RgC_mul(matrix, coordinates) : coordinates;
}

/**
 * An ideal given by a stored HNF matrix, in local coordinates.  A fractional
 * ideal is cleared of its denominator first, so that the Hermite form is
 * taken over the integers, and scaled back afterwards.
 */
static GEN
transformed_ideal(GEN matrix, GEN hnf)
{
    GEN denominator;
    GEN integral = Q_remove_denom(hnf, &denominator);
    GEN transformed = ZM_hnf(ZM_mul(matrix, integral));
    return denominator ? RgM_Rg_div(transformed, denominator) : transformed;
}

/**
 * A factored element in local coordinates.  Rational generators carry no
 * basis and are left alone.
 */
static GEN
transformed_famat(GEN matrix, GEN famat)
{
    GEN generators = gel(famat, 1);
    GEN result = cgetg(3, t_MAT);
    GEN transformed = cgetg(lg(generators), t_COL);
    for (long i = 1; i < lg(generators); ++i)
        gel(transformed, i) =
            transformed_element(matrix, gel(generators, i));
    gel(result, 1) = transformed;
    gel(result, 2) = gel(famat, 2);
    return result;
}

static GEN
compact_principal_ideal(GEN nf, GEN compact)
{
    return idealhnf0(
        nf, idealfactorback(nf, famat_idealfactor(nf, compact), NULL, 0),
        NULL);
}

static GEN
one_minus_sigma(GEN nf, GEN sigma, GEN ideal)
{
    return idealmul(
        nf, ideal, idealinv(nf, galoisapply(nf, sigma, ideal)));
}

static GEN
relative_norm_compact(GEN Lrel, GEN K, GEN compact)
{
    GEN result = gcopy(compact);
    for (long i = 1; i < lg(gel(compact, 1)); ++i)
    {
        GEN relative =
            rnfeltabstorel(Lrel, gmael(compact, 1, i));
        gmael(result, 1, i) =
            algtobasis(K, rnfeltnorm(Lrel, relative));
    }
    return result;
}

static GEN
compact_quotient(GEN numerator, GEN denominator)
{
    GEN result = cgetg(3, t_MAT);
    gel(result, 1) =
        shallowconcat(gel(numerator, 1), mkcol(denominator));
    gel(result, 2) =
        shallowconcat(gel(numerator, 2), mkcol(gen_m1));
    return result;
}

/** Validate a stored character label/vector and return its matrix slot. */
static long
character_index_and_vector(
    const char *label, GEN stored, GEN p,
    const char *entry_label, long column, GEN *normalized)
{
    static const char *labels[] = {
        "a", "b", "c", "a+b", "a+c", "b+c", "2a", "2b", "2c"
    };
    static const long coordinates[][3] = {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
        {2, 0, 0}, {0, 2, 0}, {0, 0, 2}
    };

    long index = -1;
    for (long i = 0; i < 9; ++i)
        if (strcmp(label, labels[i]) == 0)
        {
            index = i;
            break;
        }
    if (index < 0) fail(entry_label, column, "unknown character label");
    if ((typ(stored) != t_VEC && typ(stored) != t_COL)
        || glength(stored) != 3)
        fail(entry_label, column, "character vector is not length 3");

    GEN value = cgetg(4, t_COL);
    int nonzero = 0;
    for (long coordinate = 1; coordinate <= 3; ++coordinate)
    {
        if (typ(gel(stored, coordinate)) != t_INT)
            fail(entry_label, column, "character coordinate is not integral");
        gel(value, coordinate) = modii(gel(stored, coordinate), p);
        if (signe(gel(value, coordinate))) nonzero = 1;
        if (!equaliu(
                gel(value, coordinate),
                coordinates[index][coordinate - 1]))
            fail(entry_label, column, "character label/vector mismatch");
    }
    if (!nonzero) fail(entry_label, column, "character vector is zero");
    *normalized = value;
    return index;
}

/** Express PARI's relative cyclic generator in the absolute field basis. */
static GEN
relative_generator_as_absolute_automorphism(GEN Labs, GEN Lrel)
{
    GEN relative = rnfcycaut(Lrel);
    if (!relative) return NULL;
    return algtobasis(
        Labs,
        nfadd(
            Labs, rnfeltreltoabs(Lrel, relative),
            gmul(rnf_get_k(Lrel), rnf_get_alpha(Lrel))));
}

/** Find an automorphism's exponent relative to a cyclic generator. */
static long
automorphism_exponent(
    GEN Labs, GEN generator, GEN automorphism, long order)
{
    GEN absolute_generator =
        algtobasis(Labs, pol_x(nf_get_varn(Labs)));
    GEN current = absolute_generator;
    for (long exponent = 1; exponent < order; ++exponent)
    {
        current = galoisapply(Labs, generator, current);
        if (gequal(current, automorphism)) return exponent;
    }
    return 0;
}

/**
 * Compute base-class Artin exponents in both PARI's and stored sigma's
 * orientations.
 */
static GEN
artin_exponents_relative_to(
    GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN p,
    const char *label, long column, GEN *canonical_exponents,
    long *sigma_exponent)
{
    long p_int = itos(p);
    GEN relative_generator =
        relative_generator_as_absolute_automorphism(Labs, Lrel);
    if (!relative_generator)
        fail(label, column, "rnfcycaut did not return a generator");
    long sigma_power =
        automorphism_exponent(Labs, relative_generator, sigma, p_int);
    if (!sigma_power)
        fail(label, column, "stored sigma is not a power of rnfcycaut");

    GEN relative_automorphism = rnfcycaut(Lrel);
    GEN all_automorphisms = allauts(Lrel, relative_automorphism);
    GEN generators = bnf_get_gen(K);
    GEN cyc = bnf_get_cyc(K);
    GEN canonical = cgetg(4, t_COL);
    GEN normalized = cgetg(4, t_COL);
    long coordinate = 1;
    GEN inverse_sigma_power = Fp_inv(stoi(sigma_power), p);

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
        gel(canonical, coordinate) = exponent;
        gel(normalized, coordinate) =
            Fp_mul(exponent, inverse_sigma_power, p);
        ++coordinate;
    }
    if (coordinate != 4)
        fail(label, column, "base p-relevant class rank is not 3");

    *canonical_exponents = canonical;
    *sigma_exponent = sigma_power;
    return normalized;
}

/** Prove that the extension is ker(character) and sigma realizes character. */
static void
verify_character_and_normalization(
    GEN Labs, GEN Lrel, GEN K, GEN sigma, GEN character, GEN p,
    const char *label, long column)
{
    GEN canonical;
    long sigma_exponent;
    GEN normalized = artin_exponents_relative_to(
        Labs, Lrel, K, sigma, p, label, column,
        &canonical, &sigma_exponent);

    long pivot = 0;
    for (long coordinate = 1; coordinate <= 3; ++coordinate)
        if (signe(gel(character, coordinate)))
        {
            pivot = coordinate;
            break;
        }
    if (!pivot) fail(label, column, "character vector is zero");

    GEN line_scale =
        Fp_div(gel(canonical, pivot), gel(character, pivot), p);
    if (!signe(line_scale))
        fail(label, column, "extension Artin character is zero");
    for (long coordinate = 1; coordinate <= 3; ++coordinate)
    {
        GEN expected_canonical =
            Fp_mul(line_scale, gel(character, coordinate), p);
        if (!equalii(gel(canonical, coordinate), expected_canonical))
            fail(label, column, "extension Artin character is not ker(x)");
        if (!equalii(gel(normalized, coordinate), gel(character, coordinate)))
            fail(label, column, "stored sigma normalization mismatch");
    }
}

/** Check that the stored absolute automorphism fixes K and has order five. */
static void
verify_automorphism(
    GEN Labs, GEN Lrel, GEN K, GEN sigma,
    const char *label, long column)
{
    GEN base_generator =
        rnfeltup0(
            Lrel, pol_x(nf_get_varn(bnf_get_nf(K))), 1);
    if (!gequal(
            galoisapply(Labs, sigma, base_generator),
            base_generator))
        fail(label, column, "stored sigma does not fix K");

    GEN absolute_generator =
        algtobasis(Labs, pol_x(nf_get_varn(Labs)));
    GEN current = absolute_generator;
    for (long exponent = 1; exponent <= 5; ++exponent)
    {
        current = galoisapply(Labs, sigma, current);
        if (exponent < 5 && gequal(current, absolute_generator))
            fail(label, column, "stored sigma has order less than 5");
    }
    if (!gequal(current, absolute_generator))
        fail(label, column, "stored sigma does not have order 5");
}

/** Resolve the unit ambiguity in AC2 by reduction at the stored odd prime. */
static void
verify_modular_sign(
    GEN K, GEN compact, GEN ell, GEN prime,
    const char *label, long column)
{
    if (!uisprime(itou(ell)) || !mpodd(ell))
        fail(label, column, "stored modular prime is not odd and prime");

    GEN decomposition = idealprimedec(K, ell);
    int found = 0;
    for (long i = 1; i < lg(decomposition); ++i)
        if (gequal(gel(decomposition, i), prime))
        {
            found = 1;
            break;
        }
    if (!found)
        fail(label, column, "stored prime ideal is not above stored ell");

    for (long i = 1; i < lg(gel(compact, 1)); ++i)
        if (nfval(K, gmael(compact, 1, i), prime) != 0)
            fail(label, column, "compact AC2 factor is not a unit at q");

    GEN modpr = nfmodprinit(K, prime);
    GEN residue = nfmodpr(K, gen_1, modpr);
    for (long i = 1; i < lg(gel(compact, 1)); ++i)
    {
        GEN factor = nfmodpr(K, gmael(compact, 1, i), modpr);
        residue =
            FF_mul(residue, FF_pow(factor, gmael(compact, 2, i)));
    }

    GEN plus_one = nfmodpr(K, gen_1, modpr);
    GEN minus_one = nfmodpr(K, gen_m1, modpr);
    if (FF_equal(plus_one, minus_one))
        fail(label, column, "+1 and -1 coincide at stored odd q");
    if (!FF_equal(residue, plus_one))
        fail(label, column, "N(t_AC)/a' is not +1 modulo q");
}

/** Return the p-relevant coordinates of an ideal in PARI's class basis. */
static GEN
class_coordinates_mod_p(GEN K, GEN ideal, GEN p)
{
    GEN exponents = bnfisprincipal0(K, ideal, 0);
    GEN cyc = bnf_get_cyc(K);
    long rank = 0;
    for (long i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p)) ++rank;

    GEN result = cgetg(rank + 1, t_COL);
    long coordinate = 1;
    for (long i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
            gel(result, coordinate++) = modii(gel(exponents, i), p);
    return result;
}

/** Read the GP expression without asking GP to execute a source file. */
static GEN
read_certificate(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (!file) fail(NULL, 0, "cannot open certificate file");
    if (fseek(file, 0, SEEK_END) != 0)
        fail(NULL, 0, "cannot seek certificate file");
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0)
        fail(NULL, 0, "cannot size certificate file");

    char *text = malloc((size_t)size + 1);
    if (!text) fail(NULL, 0, "cannot allocate certificate buffer");
    if (fread(text, 1, (size_t)size, file) != (size_t)size)
        fail(NULL, 0, "cannot read complete certificate file");
    fclose(file);
    text[size] = '\0';

    GEN certificate = gp_read_str(text);
    free(text);
    return certificate;
}

int
main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "certificate.gp";
    pari_init_opts(
        1L << 30, 1048576,
        INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    paristack_setsize(1L << 30, 1L << 33);

    /* Validate the certificate schema and certified base-field conventions. */
    GEN certificate = read_certificate(path);
    if (typ(certificate) != t_VEC || lg(certificate) != 8)
        fail(NULL, 0, "invalid top-level certificate schema");
    if (typ(gel(certificate, CERT_FORMAT)) != t_INT)
        fail(NULL, 0, "unsupported certificate format");
    long format = itos(gel(certificate, CERT_FORMAT));
    if (format != 1 && format != 2)
        fail(NULL, 0, "unsupported certificate format");
    if (!equaliu(gel(certificate, CERT_PARI_VERSION), PARI_VERSION_CODE))
        fail(NULL, 0, "PARI version differs from certificate generator");

    GEN p = gel(certificate, CERT_P);
    if (!equaliu(p, 5))
        fail(NULL, 0, "certificate prime is not 5");

    GEN K = Buchall(
        gel(certificate, CERT_BASE_POLYNOMIAL), nf_FORCE, DEFAULTPREC);
    if (bnfcertify0(K, 0) != 1)
        fail(NULL, 0, "bnfcertify(K) failed");
    if (!gequal(
            nf_get_disc(bnf_get_nf(K)),
            gel(certificate, CERT_DISCRIMINANT)))
        fail(NULL, 0, "base discriminant mismatch");

    GEN base_data = gel(certificate, CERT_BASE_DATA);
    long base_fields = format == 1 ? 6 : 7;
    if (typ(base_data) != t_VEC || lg(base_data) != base_fields)
        fail(NULL, 0, "invalid base metadata");
    if (!gequal(bnf_get_cyc(K), gel(base_data, 1))
        || !gequal(bnf_get_no(K), gel(base_data, 2))
        || !gequal(bnf_get_gen(K), gel(base_data, 3)))
        fail(NULL, 0, "certified base class-group convention mismatch");
    long p_divisible = 0;
    GEN base_cyc = bnf_get_cyc(K);
    for (long i = 1; i < lg(base_cyc); ++i)
        if (dvdii(gel(base_cyc, i), p)) ++p_divisible;
    if (p_divisible != 3)
        fail(NULL, 0, "certified base p-class rank is not 3");
    if (glength(bnf_get_fu(K)) != 0
        || bnf_get_tuN(K) != itos(gel(base_data, 4))
        || !gequal(bnf_get_tuU(K), gel(base_data, 5))
        || bnf_get_tuN(K) != 2
        || !gequal(bnf_get_tuU(K), gen_m1))
        fail(NULL, 0, "certified base unit data mismatch");

    /*
     * Format 1 carries no integral basis, so its coordinates are only
     * meaningful on a machine whose nfinit happens to reproduce the basis of
     * the generator.  Say so rather than let a failure look like bad
     * arithmetic.
     */
    GEN base_change = NULL;
    if (format == 1)
        pari_fprintf(
            stderr,
            "WARNING: certificate format 1 does not record the integral "
            "basis it was written against.  A failure below may mean that "
            "this machine's nfinit chose a different basis, not that the "
            "arithmetic is wrong.  Convert to format 2.\n");
    else
        base_change = basis_change_matrix(
            bnf_get_nf(K), gel(base_data, 6), NULL, 0);

    pari_printf("BASE_BNF_CERTIFIED=PASS\n\n");

    GEN entries = gel(certificate, CERT_ENTRIES);
    if (typ(entries) != t_VEC
        || (glength(entries) != 27 && glength(entries) != 18))
        fail(NULL, 0, "certificate must contain 18 or 27 entries");
    int has_doubled = glength(entries) == 27;

    GEN matrices[9];
    int seen[9][3] = {{0}};
    for (long i = 0; i < 9; ++i)
        matrices[i] = zeromatcopy(3, 3);

    /*
     * Verify every arithmetic entry independently before accepting its norm
     * class as one column of a secondary-norm matrix.
     */
    for (long entry_index = 1; entry_index < lg(entries); ++entry_index)
    {
        GEN entry = gel(entries, entry_index);
        if (typ(entry) != t_VEC || lg(entry) != (format == 1 ? 14 : 15))
            fail(NULL, 0, "invalid entry schema");

        const char *label = GSTR(gel(entry, ENTRY_CHARACTER));
        long column = itos(gel(entry, ENTRY_COLUMN));
        if (column < 1 || column > 3)
            fail(label, column, "column is outside 1..3");
        GEN character;
        long matrix_index = character_index_and_vector(
            label, gel(entry, ENTRY_CHARACTER_VECTOR), p,
            label, column, &character);
        if (seen[matrix_index][column - 1])
            fail(label, column, "duplicate certificate entry");
        seen[matrix_index][column - 1] = 1;

        GEN Lrel = rnfinit(
            bnf_get_nf(K), gel(entry, ENTRY_RELATIVE_POLYNOMIAL));
        if (!gequal(
                rnf_get_polabs(Lrel),
                gel(entry, ENTRY_ABSOLUTE_POLYNOMIAL)))
            fail(
                label, column,
                "relative/absolute field models are incompatible");
        GEN Labs =
            nfinit0(gel(entry, ENTRY_ABSOLUTE_POLYNOMIAL), 0, DEFAULTPREC);
        if (nf_get_degree(Labs) / nf_get_degree(bnf_get_nf(K)) != 5)
            fail(label, column, "relative degree is not 5");
        if (!equalii(
                nf_get_disc(Labs),
                powiu(nf_get_disc(bnf_get_nf(K)), 5)))
            fail(label, column, "relative discriminant is not trivial");

        GEN sigma = gel(entry, ENTRY_SIGMA);
        GEN I_prime = gel(entry, ENTRY_I_PRIME);
        GEN a_prime = gel(entry, ENTRY_A_PRIME);
        GEN J = gel(entry, ENTRY_J);
        GEN t_AC = gel(entry, ENTRY_T_AC);

        /*
         * Translate the stored coordinates into this machine's bases before
         * anything is checked, so that every test below sees the algebraic
         * objects the generator meant.  Where the two bases agree the matrix
         * is the identity and this changes nothing.
         */
        if (format != 1)
        {
            GEN change = basis_change_matrix(
                Labs, gel(entry, ENTRY_INTEGRAL_BASIS), label, column);
            sigma = transformed_element(change, sigma);
            I_prime = transformed_ideal(change, I_prime);
            t_AC = transformed_famat(change, t_AC);
            a_prime = transformed_element(base_change, a_prime);
            J = transformed_ideal(base_change, J);
        }

        verify_automorphism(Labs, Lrel, K, sigma, label, column);
        verify_character_and_normalization(
            Labs, Lrel, K, sigma, character, p, label, column);
        GEN unit_L = idealhnf0(Labs, gen_1, NULL);
        GEN unit_K = idealhnf0(K, gen_1, NULL);

        GEN base_relation = idealmul(
            K, idealhnf0(K, a_prime, NULL), idealpow(K, J, p));
        if (!gequal(idealhnf0(K, base_relation, NULL), unit_K))
            fail(label, column, "(a') J^5 is not O_K");

        GEN operated = one_minus_sigma(Labs, sigma, I_prime);
        operated = one_minus_sigma(Labs, sigma, operated);
        GEN t_ideal = compact_principal_ideal(Labs, t_AC);
        GEN iJ = rnfidealup0(Lrel, J, 1);
        GEN ac1 = idealmul(
            Labs, idealmul(Labs, operated, t_ideal), iJ);
        if (!gequal(idealhnf0(Labs, ac1, NULL), unit_L))
            fail(
                label, column,
                "(1-sigma)^2 I' (t_AC) i(J) is not O_L");

        GEN norm_t = relative_norm_compact(Lrel, K, t_AC);
        GEN quotient = compact_quotient(norm_t, a_prime);
        if (!gequal(compact_principal_ideal(K, quotient), unit_K))
            fail(label, column, "N(t_AC)/a' does not generate O_K");
        verify_modular_sign(
            K, quotient, gel(entry, ENTRY_ELL), gel(entry, ENTRY_PRIME),
            label, column);

        GEN relative_I = rnfidealabstorel(Lrel, I_prime);
        GEN norm_I = rnfidealnormrel(Lrel, relative_I);
        GEN coordinates = class_coordinates_mod_p(K, norm_I, p);
        if (!gequal(coordinates, gel(entry, ENTRY_NORM_CLASS)))
            fail(label, column, "norm-class coordinates mismatch");
        gel(matrices[matrix_index], column) = gcopy(coordinates);

        pari_printf(
            "%s/e%ld  FIELD_MODEL_COMPATIBILITY=PASS "
            "ARTIN_CHARACTER=PASS SIGMA_NORMALIZATION=PASS "
            "AC1=PASS AC2=PASS NORM_CLASS=%Ps\n",
            label, column, gtovec(coordinates));
    }

    static const char *all_labels[] = {
        "a", "b", "c", "a+b", "a+c", "b+c", "2a", "2b", "2c"
    };
    for (long character = 0; character < (has_doubled ? 9 : 6); ++character)
        for (long column = 0; column < 3; ++column)
            if (!seen[character][column])
                fail(
                    all_labels[character], column + 1,
                    "missing certificate entry");

    /*
     * For the principal example the six matrices are compared with the
     * published values; for any other field they are printed, and may be
     * compared against a result record passed as the second argument.
     */
    static const char *matrix_labels[] = {
        "D_a", "D_b", "D_c", "D_(a+b)", "D_(a+c)", "D_(b+c)"
    };
    if (equalii(
            gel(certificate, CERT_DISCRIMINANT), stoi(-11203620)))
    {
        static const long expected[6][3][3] = {
            {{0,0,0},{0,3,3},{0,1,0}},
            {{0,0,4},{1,0,2},{2,0,4}},
            {{1,4,0},{2,0,0},{0,1,0}},
            {{0,0,0},{3,2,1},{1,4,2}},
            {{1,1,4},{4,4,1},{0,4,0}},
            {{3,0,0},{1,3,2},{2,2,3}}
        };
        for (long matrix = 0; matrix < 6; ++matrix)
        {
            for (long column = 1; column <= 3; ++column)
                for (long row = 1; row <= 3; ++row)
                    if (!equaliu(
                            gcoeff(matrices[matrix], row, column),
                            expected[matrix][column - 1][row - 1]))
                        fail(NULL, 0, "reconstructed main matrix mismatch");
            pari_printf("%s=PASS\n", matrix_labels[matrix]);
        }
    }
    else
        for (long matrix = 0; matrix < 6; ++matrix)
            pari_printf(
                "%s=%Ps\n", matrix_labels[matrix], matrices[matrix]);

    if (argc > 2)
    {
        GEN record = read_certificate(argv[2]);
        GEN samples = NULL;
        if (typ(record) == t_VEC)
            for (long i = 1; i < lg(record); ++i)
            {
                GEN pair = gel(record, i);
                if (typ(pair) == t_VEC && lg(pair) == 3
                    && typ(gel(pair, 1)) == t_STR
                    && strcmp(
                        GSTR(gel(pair, 1)),
                        "secondary_norm_samples") == 0)
                    samples = gel(pair, 2);
            }
        if (samples == NULL || typ(samples) != t_VEC
            || glength(samples) != 6)
            fail(NULL, 0, "result record has no secondary-norm samples");
        for (long matrix = 0; matrix < 6; ++matrix)
            for (long column = 1; column <= 3; ++column)
                for (long row = 1; row <= 3; ++row)
                    if (!equalii(
                            gcoeff(matrices[matrix], row, column),
                            modii(
                                gcoeff(
                                    gel(samples, matrix + 1),
                                    row, column), p)))
                        fail(
                            NULL, 0,
                            "certificate disagrees with result record");
        pari_printf("RESULT_RECORD_MATCH=PASS\n");
    }

    /* Check the three redundant doubled-character homogeneity witnesses. */
    if (has_doubled)
        for (long matrix = 0; matrix < 3; ++matrix)
        {
            for (long column = 1; column <= 3; ++column)
                for (long row = 1; row <= 3; ++row)
                    if (!equalii(
                            gcoeff(matrices[6 + matrix], row, column),
                            Fp_mul(
                                stoi(4),
                                gcoeff(matrices[matrix], row, column), p)))
                        fail(NULL, 0, "doubled-character scaling mismatch");
            pari_printf(
                "D_(2%s)=4D_%s PASS\n",
                matrix == 0 ? "a" : matrix == 1 ? "b" : "c",
                matrix == 0 ? "a" : matrix == 1 ? "b" : "c");
        }

    pari_printf("\nCERTIFICATE VERIFIED\n");
    pari_close();
    return EXIT_SUCCESS;
}
