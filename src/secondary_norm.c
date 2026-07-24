// MIT License

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pari/pari.h>
#include "../headers/secondary_norm.h"
#include "../headers/artin_symbol.h"
#include "../headers/ext_and_aut.h"
#include "../headers/misc_functions.h"
#include "../headers/pari_internal.h"

static int
secondary_norm_diagnostics_enabled(void)
{
    const char *value = getenv("MASSEY_DIAGNOSTICS");
    const char *rank3 = getenv("MASSEY_RANK3_TEST");
    const char *audit = getenv("MASSEY_ARITHMETIC_AUDIT");
    return (value && strcmp(value, "1") == 0)
        || (rank3 && strcmp(rank3, "1") == 0)
        || (audit && strcmp(audit, "1") == 0);
}

static void
secondary_norm_error(const char *message)
{
    pari_err(e_MISC, "my_secondary_norm_operator: %s", message);
}

long
my_p_rank_from_cyc(GEN cyc, GEN p)
{
    long i, rank = 0;

    if (typ(cyc) != t_VEC && typ(cyc) != t_COL)
        pari_err_TYPE("my_p_rank_from_cyc [cyclic invariants]", cyc);
    if (typ(p) != t_INT)
        pari_err_TYPE("my_p_rank_from_cyc [prime]", p);

    for (i = 1; i < lg(cyc); ++i)
    {
        if (typ(gel(cyc, i)) != t_INT)
            pari_err_TYPE(
                "my_p_rank_from_cyc [cyclic invariant]", gel(cyc, i));
        if (dvdii(gel(cyc, i), p)) ++rank;
    }
    return rank;
}

long
my_p_class_rank(GEN K, GEN p)
{
    return my_p_rank_from_cyc(bnf_get_cyc(K), p);
}

GEN
my_p_relevant_coordinates(GEN K, GEN values, GEN p)
{
    GEN cyc = bnf_get_cyc(K);
    long i, j = 1, p_rk = my_p_class_rank(K, p);

    GEN result = cgetg(p_rk + 1, t_COL);
    for (i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
            gel(result, j++) = modii(gel(values, i), p);
    return result;
}

GEN
my_subgroup_character(GEN K, GEN H, GEN p)
{
    pari_sp av = avma;
    GEN cyc = bnf_get_cyc(K);
    long i, j, row, p_rk = my_p_class_rank(K, p);

    GEN H_p = cgetg(lg(H), t_MAT);
    for (j = 1; j < lg(H); ++j)
    {
        GEN column = cgetg(p_rk + 1, t_COL);
        row = 1;
        for (i = 1; i < lg(cyc); ++i)
            if (dvdii(gel(cyc, i), p))
                gel(column, row++) = modii(gcoeff(H, i, j), p);
        gel(H_p, j) = column;
    }

    GEN kernel = FpM_ker(shallowtrans(H_p), p);
    if (lg(kernel) != 2)
        secondary_norm_error(
            "subgroup character nullspace is not one-dimensional");

    GEN c = gcopy(gel(kernel, 1));
    GEN first = NULL;
    for (i = 1; i < lg(c); ++i)
        if (signe(gel(c, i))) { first = gel(c, i); break; }
    if (!first)
        secondary_norm_error("subgroup character is zero");

    GEN inverse = Fp_inv(first, p);
    for (i = 1; i < lg(c); ++i)
        gel(c, i) = Fp_mul(gel(c, i), inverse, p);
    return gerepilecopy(av, c);
}

long
my_sigma_exponent(GEN Labs, GEN Lrel, GEN sigma_H90, GEN p)
{
    pari_sp av = avma;
    GEN sigma_rel = rnfcycaut(Lrel);
    if (!sigma_rel)
        secondary_norm_error("rnfcycaut did not return an Artin generator");

    GEN sigma_abs_alg =
        nfadd(Labs, rnfeltreltoabs(Lrel, sigma_rel),
              gmul(rnf_get_k(Lrel), rnf_get_alpha(Lrel)));
    GEN sigma_abs = algtobasis(Labs, sigma_abs_alg);
    GEN abs_generator =
        algtobasis(Labs, pol_x(nf_get_varn(bnf_get_nf(Labs))));
    GEN current = abs_generator;
    long a, p_int = itos(p);

    for (a = 1; a < p_int; ++a)
    {
        current = galoisapply(Labs, sigma_abs, current);
        if (gequal(current, sigma_H90))
            return gc_long(av, a);
    }
    secondary_norm_error(
        "H90 automorphism is not a nonzero power of the Artin automorphism");
    return 0;
}

GEN
my_automorphism_power_checked(
    GEN Labs, GEN Lrel, GEN K, GEN sigma, long exponent, GEN p)
{
    pari_sp av = avma;
    GEN abs_generator =
        algtobasis(Labs, pol_x(nf_get_varn(bnf_get_nf(Labs))));
    GEN sigma_power = abs_generator;
    long i, p_int = itos(p);

    if (exponent <= 0 || exponent >= p_int)
        secondary_norm_error("automorphism exponent must lie in F_p^*");

    for (i = 1; i <= exponent; ++i)
        sigma_power = galoisapply(Labs, sigma, sigma_power);

    if (exponent == 2)
    {
        GEN composed = galoisapply(Labs, sigma, sigma);
        if (!gequal(sigma_power, composed))
            secondary_norm_error(
                "sigma squared does not equal sigma composed with sigma");
    }

    GEN base_generator =
        rnfeltup0(Lrel, pol_x(nf_get_varn(bnf_get_nf(K))), 1);
    if (!gequal(galoisapply(Labs, sigma_power, base_generator),
                base_generator))
        secondary_norm_error("automorphism power does not fix the base field");
    if (gequal(sigma_power, abs_generator))
        secondary_norm_error("automorphism power is trivial");

    GEN current = abs_generator;
    for (i = 1; i <= p_int; ++i)
    {
        current = galoisapply(Labs, sigma_power, current);
        if (i < p_int && gequal(current, abs_generator))
            secondary_norm_error(
                "automorphism power does not generate Gal(L/K)");
    }
    if (!gequal(current, abs_generator))
        secondary_norm_error("automorphism power does not have order p");

    return gerepilecopy(av, sigma_power);
}

static int
secondary_norm_arithmetic_audit_enabled(void)
{
    const char *value = getenv("MASSEY_ARITHMETIC_AUDIT");
    return value && strcmp(value, "1") == 0;
}

static GEN
secondary_norm_compact_inverse(GEN compact)
{
    GEN inverse = cgetg(3, t_MAT);
    gel(inverse, 1) = gcopy(gel(compact, 1));
    gel(inverse, 2) = gneg(gel(compact, 2));
    return inverse;
}

/*
 * Exact principal fractional ideal of a compact nf_GENMAT element.  This
 * never expands the algebraic element: each factor is converted to its
 * principal ideal and raised to its signed exponent.
 */
static GEN
secondary_norm_compact_principal_ideal(GEN nf, GEN compact)
{
    pari_sp av = avma;
    /*
     * PARI 2.17.4 documents famat_idealfactor as the coefficient-explosion
     * free way to obtain the ideal generated by a compact element.  It
     * factors the principal ideal of every component and combines the signed
     * exponents before idealfactorback materializes the final HNF.
     */
    GEN factorization = famat_idealfactor(nf, compact);
    GEN result = idealfactorback(nf, factorization, NULL, 0);
    return gerepilecopy(av, idealhnf0(nf, result, NULL));
}

/*
 * Evaluate a compact K-element exactly in a residue field.  We choose an odd
 * prime ideal at which every individual factor has valuation zero, so all
 * signed powers (and all represented denominators) reduce and are invertible.
 * Return [rational prime, prime ideal, residue].
 */
static GEN
secondary_norm_compact_modular_sign(GEN K, GEN compact)
{
    pari_sp av = avma;
    for (ulong ell = 3; ; ell += 2)
    {
        if (!uisprime(ell)) continue;
        GEN ell_gen = utoipos(ell);
        GEN primes = idealprimedec(K, ell_gen);
        for (long q = 1; q < lg(primes); ++q)
        {
            GEN prime = gel(primes, q);
            int suitable = 1;
            for (long i = 1; i < lg(gel(compact, 1)); ++i)
                if (nfval(K, gmael(compact, 1, i), prime) != 0)
                {
                    suitable = 0;
                    break;
                }
            if (!suitable) continue;

            GEN modpr = nfmodprinit(K, prime);
            GEN residue = nfmodpr(K, gen_1, modpr);
            for (long i = 1; i < lg(gel(compact, 1)); ++i)
            {
                GEN factor =
                    nfmodpr(K, gmael(compact, 1, i), modpr);
                GEN factor_power =
                    FF_pow(factor, gmael(compact, 2, i));
                residue = FF_mul(residue, factor_power);
            }
            return gerepilecopy(
                av, mkvec3(ell_gen, prime, residue));
        }
    }
}

static GEN
secondary_norm_audit_column(
    GEN Labs, GEN Lrel, GEN K, GEN sigma_t, GEN p,
    GEN Ja, GEN I_prime, GEN production_column,
    long input_index)
{
    pari_sp av = avma;
    GEN a_prime = gel(Ja, 1);
    GEN J = gel(Ja, 2);
    GEN iJ = rnfidealup0(Lrel, J, 1);
    GEN operated = gcopy(I_prime);
    operated = my_1MS_ideal(Labs, sigma_t, operated);
    operated = my_1MS_ideal(Labs, sigma_t, operated);

    /*
     * For n=2 the production solver makes iJ * (1-sigma)^2 I principal.
     * Its compact generator t_code is inverse to the t in AC1/AC2.
     */
    GEN principal_ideal = idealmul(Labs, iJ, operated);
    GEN principal_data =
        bnfisprincipal0(Labs, principal_ideal, nf_GENMAT);
    if (!ZV_equal0(gel(principal_data, 1)))
        secondary_norm_error("audit AC1 ideal is not principal");
    GEN t_code = gel(principal_data, 2);

    GEN norm_code = my_rel_norm_compact(Labs, Lrel, K, t_code);
    GEN norm_times_a = cgetg(3, t_MAT);
    gel(norm_times_a, 1) =
        shallowconcat(gel(norm_code, 1), mkcol(a_prime));
    gel(norm_times_a, 2) =
        shallowconcat(gel(norm_code, 2), mkcol(gen_1));
    GEN unit_exp = bnfisunit0(K, norm_times_a, NULL);
    if (glength(unit_exp) == 0)
        secondary_norm_error("audit N(t_code)*a' is not a base unit");

    GEN corrected_code = gcopy(t_code);
    if (!ZV_equal0(unit_exp))
    {
        GEN norm_operator = my_norm_operator(Labs, Lrel, K, p);
        GEN unit_solution =
            matsolvemod(
                norm_operator, zerocol(glength(unit_exp)),
                gtocol(unit_exp), 0);
        if (gequal0(unit_solution))
            secondary_norm_error(
                "audit base-unit discrepancy is not an extension-unit norm");
        GEN extension_units =
            shallowconcat(bnf_get_fu(Labs), bnf_get_tuU(Labs));
        GEN correction = cgetg(3, t_MAT);
        gel(correction, 1) = gtocol(extension_units);
        gel(correction, 2) = gneg(unit_solution);
        corrected_code = concatenate_rows(t_code, correction);
    }

    GEN t_AC = secondary_norm_compact_inverse(corrected_code);
    GEN norm_AC = my_rel_norm_compact(Labs, Lrel, K, t_AC);
    GEN norm_AC_over_a = cgetg(3, t_MAT);
    gel(norm_AC_over_a, 1) =
        shallowconcat(gel(norm_AC, 1), mkcol(a_prime));
    gel(norm_AC_over_a, 2) =
        shallowconcat(gel(norm_AC, 2), mkcol(gen_m1));
    GEN t_AC_ideal =
        secondary_norm_compact_principal_ideal(Labs, t_AC);
    GEN ac1_product =
        idealmul(
            Labs,
            idealmul(Labs, operated, t_AC_ideal),
            iJ);
    GEN ac1_hnf = idealhnf0(Labs, ac1_product, NULL);
    GEN unit_hnf = idealhnf0(Labs, gen_1, NULL);
    if (!gequal(ac1_hnf, unit_hnf))
        secondary_norm_error(
            "audit AC1 compact-element divisor equality failed");

    GEN base_relation =
        idealmul(
            K, idealhnf0(K, a_prime, NULL),
            idealpow(K, J, p));
    GEN base_unit_hnf = idealhnf0(K, gen_1, NULL);
    if (!gequal(idealhnf0(K, base_relation, NULL), base_unit_hnf))
        secondary_norm_error("audit div(a')+pJ identity failed");

    /*
     * AC1_COMPACT_DIVISOR plus div(a')+pJ=0 proves that the compact quotient
     * u=N(t_AC)/a' has zero divisor.  Confirm the resulting principal ideal
     * directly as an additional exact ideal-arithmetic check.
     */
    GEN quotient_ideal =
        secondary_norm_compact_principal_ideal(K, norm_AC_over_a);
    if (!gequal(quotient_ideal, base_unit_hnf))
        secondary_norm_error("audit AC2 quotient does not have zero divisor");

    GEN modular = secondary_norm_compact_modular_sign(K, norm_AC_over_a);
    GEN ell = gel(modular, 1);
    GEN prime = gel(modular, 2);
    GEN residue = gel(modular, 3);
    GEN modpr = nfmodprinit(K, prime);
    GEN residue_one = nfmodpr(K, gen_1, modpr);
    GEN residue_minus_one = nfmodpr(K, gen_m1, modpr);
    if (FF_equal(residue_one, residue_minus_one))
        secondary_norm_error("audit odd-prime residues do not distinguish signs");
    if (!FF_equal(residue, residue_one))
        secondary_norm_error("audit AC2 modular sign is not +1");

    /*
     * Secondary cross-check only.  For this imaginary quadratic K the sole
     * coordinate is the exponent of the torsion generator -1 modulo 2.
     */
    GEN ac2_unit_exp = bnfisunit0(K, norm_AC_over_a, NULL);
    if (glength(ac2_unit_exp) == 0 || !ZV_equal0(ac2_unit_exp))
        secondary_norm_error("audit AC2 unit coordinate is not zero");

    GEN I_rel = rnfidealabstorel(Lrel, I_prime);
    GEN norm_ideal = rnfidealnormrel(Lrel, I_rel);
    GEN class_exp = bnfisprincipal0(K, norm_ideal, 0);
    GEN cyc = bnf_get_cyc(K);
    long rank = my_p_class_rank(K, p), coordinate = 1;
    GEN independent = cgetg(rank + 1, t_COL);
    for (long i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
            gel(independent, coordinate++) = modii(gel(class_exp, i), p);
    if (!gequal(independent, production_column))
        secondary_norm_error(
            "audit independent norm-class extraction disagrees with production");

    pari_printf(
        "  input e_%ld: AC1_COMPACT_DIVISOR=PASS "
        "AC2_UNIT_STATUS=PASS AC2_MODULAR_SIGN=PASS "
        "AC2_UNIT_COORD=PASS "
        "independent D=%Ps production D=%Ps MATCH\n",
        input_index, gtovec(independent), gtovec(production_column));
    pari_printf(
        "    rational_prime=%Ps prime_ideal=%Ps residue(u)=%Ps "
        "residue(+1)=%Ps residue(-1)=%Ps unit_coordinates=%Ps\n",
        ell, prime, residue, residue_one, residue_minus_one,
        ac2_unit_exp);
    return gerepilecopy(
        av, mkvec4(t_AC, I_prime, norm_ideal, independent));
}

GEN
my_secondary_norm_operator(
    GEN K, GEN p, GEN prescribed_character_t,
    GEN Ja_vect, GEN D_prime_vect)
{
    pari_sp av = avma;
    GEN cyc = bnf_get_cyc(K);
    long i, j, p_int = itos(p), p_rk = my_p_class_rank(K, p);

    if (p_int <= 3)
        secondary_norm_error("p <= 3 is not supported");
    if (typ(prescribed_character_t) != t_VEC
        && typ(prescribed_character_t) != t_COL)
        pari_err_TYPE(
            "my_secondary_norm_operator [prescribed character]",
            prescribed_character_t);

    if (glength(prescribed_character_t) != p_rk)
        pari_err_DIM("my_secondary_norm_operator [prescribed character]");

    GEN t = cgetg(p_rk + 1, t_COL);
    int nonzero = 0;
    for (i = 1; i <= p_rk; ++i)
    {
        gel(t, i) = modii(gel(prescribed_character_t, i), p);
        if (signe(gel(t, i))) nonzero = 1;
    }
    if (!nonzero)
        secondary_norm_error("prescribed character must be nonzero");

    GEN chi = cgetg(lg(cyc), t_VEC);
    j = 1;
    for (i = 1; i < lg(cyc); ++i)
    {
        if (dvdii(gel(cyc, i), p))
        {
            gel(chi, i) =
                modii(mulii(diviiexact(gel(cyc, i), p), gel(t, j)), gel(cyc, i));
            ++j;
        }
        else gel(chi, i) = gen_0;
    }

    GEN H = charker(cyc, chi);
    GEN line = my_subgroup_character(K, H, p);
    long pivot = 0;
    for (i = 1; i <= p_rk; ++i)
        if (signe(gel(t, i))) { pivot = i; break; }
    GEN line_scale = Fp_div(gel(line, pivot), gel(t, pivot), p);
    for (i = 1; i <= p_rk; ++i)
        if (!equalii(gel(line, i), Fp_mul(line_scale, gel(t, i), p)))
            secondary_norm_error(
                "charker subgroup annihilator is not the prescribed line");

    GEN classfield_polynomials =
        bnrclassfield(K, H, 0, DEFAULTPREC);
    GEN K_ext =
        my_ext(K, mkvec(classfield_polynomials), p, 1, D_prime_vect);
    GEN Labs = gmael(K_ext, 1, 1);
    GEN Lrel = gmael(K_ext, 1, 2);
    GEN sigma_H90 = gmael(K_ext, 1, 3);
    GEN Lbnr = gmael(K_ext, 1, 4);

    (void)rnfidealup0(Lrel, idealhnf0(K, gen_1, NULL), 1);

    GEN A = zerocol(p_rk);
    GEN generators = bnf_get_gen(K);
    j = 1;
    for (i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
        {
            long value =
                my_Artin_symbol(
                    Labs, Lrel, K, gel(generators, i), p_int);
            gel(A, j++) = stoi((value % p_int + p_int) % p_int);
        }

    GEN lambda = Fp_div(gel(A, pivot), gel(t, pivot), p);
    if (!signe(lambda))
        secondary_norm_error("Artin character is zero");
    for (i = 1; i <= p_rk; ++i)
        if (!equalii(gel(A, i), Fp_mul(lambda, gel(t, i), p)))
            secondary_norm_error(
                "Artin character is not a scalar multiple of prescribed t");

    GEN u = stoi(my_sigma_exponent(Labs, Lrel, sigma_H90, p));
    GEN a = Fp_div(lambda, u, p);
    GEN sigma_t =
        my_automorphism_power_checked(
            Labs, Lrel, K, sigma_H90, itos(a), p);
    GEN sigma_t_exponent =
        stoi(my_sigma_exponent(Labs, Lrel, sigma_t, p));
    if (!equalii(sigma_t_exponent, lambda))
        secondary_norm_error(
            "normalized automorphism exponent does not equal lambda");
    for (i = 1; i <= p_rk; ++i)
        if (!equalii(
                Fp_mul(Fp_inv(sigma_t_exponent, p), gel(A, i), p),
                gel(t, i)))
            secondary_norm_error(
                "normalized H90 automorphism does not represent prescribed t");

    GEN I_prime_vect =
        my_H90_vect_2(
            Labs, Lrel, Lbnr, K, sigma_t, Ja_vect, p, 2);
    long inputs = glength(Ja_vect);
    GEN D = cgetg(inputs + 1, t_MAT);
    for (j = 1; j <= inputs; ++j)
    {
        GEN I_rel =
            rnfidealabstorel(Lrel, gel(I_prime_vect, j));
        GEN NI = rnfidealnormrel(Lrel, I_rel);
        GEN class_exp = bnfisprincipal0(K, NI, 0);
        gel(D, j) = my_p_relevant_coordinates(K, class_exp, p);
    }

    if (secondary_norm_diagnostics_enabled())
    {
        pari_printf("\nMASSEY_DIAGNOSTICS prescribed secondary norm\n");
        pari_printf("  prescribed t = %Ps\n", gtovec(t));
        pari_printf("  PARI character chi = %Ps\n", chi);
        pari_printf("  subgroup H = %Ps\n", H);
        pari_printf("  Artin character A = %Ps\n", gtovec(A));
        pari_printf("  lambda = %Ps\n", lambda);
        pari_printf("  sigma_H90 exponent u = %Ps\n", u);
        pari_printf("  normalization exponent a = lambda/u = %Ps\n", a);
        pari_printf(
            "  verified sigma_t exponent relative to sigma_Artin = %Ps\n",
            sigma_t_exponent);
        pari_printf("  D_t = %Ps\n", D);
    }

    if (secondary_norm_arithmetic_audit_enabled())
    {
        GEN Hdet = det(H);
        GEN Kdisc = nf_get_disc(bnf_get_nf(K));
        GEN Ldisc = nf_get_disc(bnf_get_nf(Labs));
        long relative_degree =
            nf_get_degree(bnf_get_nf(Labs)) / nf_get_degree(bnf_get_nf(K));

        if (!equalii(absi_shallow(Hdet), p))
            secondary_norm_error("audit subgroup index is not p");
        if (relative_degree != p_int)
            secondary_norm_error("audit class field does not have degree p");
        /*
         * H is passed to bnrclassfield for the ordinary class group (modulus
         * one), hence the class-field conductor has no finite prime.  The
         * discriminant identity is an independent arithmetic check that the
         * relative discriminant ideal is one.
         */
        if (!equalii(Ldisc, powiu(Kdisc, p_int)))
            secondary_norm_error("audit extension has nontrivial relative discriminant");

        pari_printf("\nMASSEY_ARITHMETIC_AUDIT character certificate\n");
        pari_printf("  prescribed t = %Ps\n", gtovec(t));
        pari_printf("  chi = %Ps\n", chi);
        pari_printf("  H = %Ps\n", H);
        pari_printf("  det(H) = %Ps; [Cl(K):H] = %Ps\n", Hdet, absi_shallow(Hdet));
        pari_printf("  class-field modulus = 1 (ordinary class group)\n");
        pari_printf(
            "  relative degree = %ld; disc(L) = disc(K)^5 = %Ps; "
            "finite-unramified=PASS\n",
            relative_degree, Ldisc);
        for (i = 1; i < lg(H); ++i)
        {
            GEN dot = gen_0;
            long coordinate = 1;
            for (j = 1; j < lg(cyc); ++j)
                if (dvdii(gel(cyc, j), p))
                    dot = Fp_add(
                        dot,
                        Fp_mul(gel(t, coordinate++), gcoeff(H, j, i), p),
                        p);
            if (signe(dot))
                secondary_norm_error("audit H column is not in ker(t)");
            pari_printf("  H column %ld: t(h)=0 PASS\n", i);
        }
        pari_printf(
            "  A=%Ps lambda=%Ps u=%Ps lambda/u=%Ps "
            "sigma_t/sigma_Artin=%Ps\n",
            gtovec(A), lambda, u, a, sigma_t_exponent);
        for (i = 1; i <= p_rk; ++i)
        {
            GEN verified =
                Fp_mul(Fp_inv(sigma_t_exponent, p), gel(A, i), p);
            pari_printf(
                "  generator g_%ld: Artin exponent A_i=%Ps; "
                "relative to sigma_t=%Ps; prescribed t_i=%Ps PASS\n",
                i, gel(A, i), verified, gel(t, i));
        }

        pari_printf("  direct AC and independent norm-class checks\n");
        for (j = 1; j <= inputs; ++j)
            (void)secondary_norm_audit_column(
                Labs, Lrel, K, sigma_t, p,
                gel(Ja_vect, j), gel(I_prime_vect, j), gel(D, j), j);
    }

    return gerepilecopy(av, D);
}

GEN
my_secondary_norm_basis_family(
    GEN K, GEN p, GEN Ja_vect, GEN D_prime_vect)
{
    pari_sp av = avma;
    long i, j, pair = 1;

    /*
     * Character-space dimension comes from Cl(K)/p.  D_prime_vect remains
     * independent discriminant-prime data used only by my_ext.
     */
    long m = my_p_class_rank(K, p);
    if (!m)
        secondary_norm_error("Cl(K)/p is trivial");

    GEN basis = cgetg(m + 1, t_VEC);
    for (i = 1; i <= m; ++i)
    {
        struct timespec started, finished;
        GEN t = zerocol(m);
        gel(t, i) = gen_1;
        clock_gettime(CLOCK_MONOTONIC, &started);
        gel(basis, i) =
            my_secondary_norm_operator(
                K, p, t, Ja_vect, D_prime_vect);
        clock_gettime(CLOCK_MONOTONIC, &finished);
        const char *rank3 = getenv("MASSEY_RANK3_TEST");
        if (rank3 && strcmp(rank3, "1") == 0)
        {
            double elapsed =
                (double)(finished.tv_sec - started.tv_sec)
                + (double)(finished.tv_nsec - started.tv_nsec) / 1e9;
            pari_printf(
                "MASSEY_RANK3_TEST sample t=%Ps elapsed=%.3f s\n",
                gtovec(t), elapsed);
            fflush(stdout);
        }
    }

    GEN pairs = cgetg((m * (m - 1)) / 2 + 1, t_VEC);
    for (i = 1; i <= m; ++i)
        for (j = i + 1; j <= m; ++j)
        {
            struct timespec started, finished;
            GEN t = zerocol(m);
            gel(t, i) = gen_1;
            gel(t, j) = gen_1;
            clock_gettime(CLOCK_MONOTONIC, &started);
            gel(pairs, pair++) =
                my_secondary_norm_operator(
                    K, p, t, Ja_vect, D_prime_vect);
            clock_gettime(CLOCK_MONOTONIC, &finished);
            const char *rank3 = getenv("MASSEY_RANK3_TEST");
            if (rank3 && strcmp(rank3, "1") == 0)
            {
                double elapsed =
                    (double)(finished.tv_sec - started.tv_sec)
                    + (double)(finished.tv_nsec - started.tv_nsec) / 1e9;
                pari_printf(
                    "MASSEY_RANK3_TEST sample t=%Ps elapsed=%.3f s\n",
                    gtovec(t), elapsed);
                fflush(stdout);
            }
        }

    GEN D0 = gel(basis, 1);
    long columns = lg(D0);
    long rows = lg(gel(D0, 1));
    for (i = 1; i <= m; ++i)
        if (typ(gel(basis, i)) != t_MAT
            || lg(gel(basis, i)) != columns
            || lg(gmael(basis, i, 1)) != rows)
            secondary_norm_error(
                "basis-family matrices have inconsistent dimensions");
    for (i = 1; i <= glength(pairs); ++i)
        if (typ(gel(pairs, i)) != t_MAT
            || lg(gel(pairs, i)) != columns
            || lg(gmael(pairs, i, 1)) != rows)
            secondary_norm_error(
                "pair-family matrices have inconsistent dimensions");

    return gerepilecopy(av, mkvec3(stoi(m), basis, pairs));
}
