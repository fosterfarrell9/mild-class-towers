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
    return (value && strcmp(value, "1") == 0)
        || (rank3 && strcmp(rank3, "1") == 0);
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
