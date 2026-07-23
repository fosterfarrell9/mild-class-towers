// MIT License

#include <pari/pari.h>
#include "../headers/massey_tensor.h"

static void
tensor_error(const char *message)
{
    pari_err(e_MISC, "massey_tensor: %s", message);
}

static long
pair_index(long m, long i, long j)
{
    if (i >= j || i < 1 || j > m)
        tensor_error("invalid pair index");
    return ((i - 1) * (2 * m - i)) / 2 + (j - i);
}

static long
word_index(long m, long i, long j, long k)
{
    if (i < 1 || i > m || j < 1 || j > m || k < 1 || k > m)
        tensor_error("invalid word index");
    return ((i - 1) * m + (j - 1)) * m + k;
}

static void
validate_matrix_dimensions(GEN D, long rows, long columns)
{
    if (typ(D) != t_MAT || lg(D) != columns + 1)
        tensor_error("secondary norm matrices have inconsistent columns");
    if (columns && lg(gel(D, 1)) != rows + 1)
        tensor_error("secondary norm matrices have inconsistent rows");
}

static void
family_dimensions(
    GEN family, long *m, long *rows, long *columns)
{
    if (typ(family) != t_VEC || lg(family) != 4)
        tensor_error("quadratic family must be [m,D_basis,D_pairs]");

    *m = itos(gel(family, 1));
    GEN basis = gel(family, 2), pairs = gel(family, 3);
    if (*m <= 0 || typ(basis) != t_VEC || glength(basis) != *m)
        tensor_error("invalid D_basis vector");
    if (typ(pairs) != t_VEC
        || glength(pairs) != (*m * (*m - 1)) / 2)
        tensor_error("invalid D_pairs vector");

    GEN D0 = gel(basis, 1);
    if (typ(D0) != t_MAT || lg(D0) <= 1)
        tensor_error("secondary norm matrix must have at least one column");
    *columns = lg(D0) - 1;
    *rows = lg(gel(D0, 1)) - 1;
    if (*rows != *m)
        tensor_error("secondary norm row count does not equal m");

    for (long i = 1; i <= *m; ++i)
        validate_matrix_dimensions(gel(basis, i), *rows, *columns);
    for (long i = 1; i <= glength(pairs); ++i)
        validate_matrix_dimensions(gel(pairs, i), *rows, *columns);
}

GEN
my_reconstruct_secondary_norm(
    GEN p, GEN quadratic_family, GEN t_input)
{
    pari_sp av = avma;
    long m, rows, columns;
    family_dimensions(quadratic_family, &m, &rows, &columns);
    if ((typ(t_input) != t_VEC && typ(t_input) != t_COL)
        || glength(t_input) != m)
        pari_err_DIM("my_reconstruct_secondary_norm [character]");

    GEN basis = gel(quadratic_family, 2);
    GEN pairs = gel(quadratic_family, 3);
    GEN t = cgetg(m + 1, t_COL);
    for (long i = 1; i <= m; ++i)
        gel(t, i) = modii(gel(t_input, i), p);

    GEN result = FpM_Fp_mul(gel(basis, 1), gen_0, p);
    for (long i = 1; i <= m; ++i)
    {
        GEN coefficient = Fp_sqr(gel(t, i), p);
        if (signe(coefficient))
            result = FpM_add(
                result,
                FpM_Fp_mul(gel(basis, i), coefficient, p),
                p);
    }

    for (long i = 1; i <= m; ++i)
        for (long j = i + 1; j <= m; ++j)
        {
            GEN coefficient = Fp_mul(gel(t, i), gel(t, j), p);
            if (signe(coefficient))
            {
                GEN cross =
                    FpM_sub(
                        FpM_sub(
                            gel(pairs, pair_index(m, i, j)),
                            gel(basis, i), p),
                        gel(basis, j), p);
                result = FpM_add(
                    result,
                    FpM_Fp_mul(cross, coefficient, p),
                    p);
            }
        }
    return gerepilecopy(av, result);
}

GEN
my_secondary_norm_delta_basis(
    GEN p, GEN quadratic_family, long i, long k)
{
    pari_sp av = avma;
    long m, rows, columns;
    family_dimensions(quadratic_family, &m, &rows, &columns);
    if (i < 1 || i > m || k < 1 || k > m)
        tensor_error("DeltaD basis index out of range");

    GEN basis = gel(quadratic_family, 2);
    GEN result;
    if (i == k)
        result =
            FpM_Fp_mul(gel(basis, i), subii(p, stoi(2)), p);
    else
    {
        long first = i < k ? i : k;
        long second = i < k ? k : i;
        GEN pair =
            gel(gel(quadratic_family, 3),
                pair_index(m, first, second));
        result =
            FpM_sub(
                FpM_add(gel(basis, i), gel(basis, k), p),
                pair, p);
    }
    return gerepilecopy(av, result);
}

GEN
my_triple_massey_word_matrix(
    GEN p, GEN quadratic_family)
{
    pari_sp av = avma;
    long m, rows, inputs;
    family_dimensions(quadratic_family, &m, &rows, &inputs);
    GEN T = cgetg(m * m * m + 1, t_MAT);

    for (long i = 1; i <= m; ++i)
        for (long j = 1; j <= m; ++j)
            for (long k = 1; k <= m; ++k)
            {
                GEN delta =
                    my_secondary_norm_delta_basis(
                        p, quadratic_family, i, k);
                GEN evaluation = cgetg(inputs + 1, t_COL);
                for (long l = 1; l <= inputs; ++l)
                    gel(evaluation, l) = gcoeff(delta, j, l);
                gel(T, word_index(m, i, j, k)) = evaluation;
            }
    return gerepilecopy(av, T);
}

GEN
my_triple_massey_contract(
    GEN p, GEN T, long m, GEN x, GEN y, GEN z)
{
    pari_sp av = avma;
    if (typ(T) != t_MAT || lg(T) != m * m * m + 1)
        tensor_error("word matrix has the wrong number of columns");
    if ((typ(x) != t_VEC && typ(x) != t_COL) || glength(x) != m
        || (typ(y) != t_VEC && typ(y) != t_COL) || glength(y) != m
        || (typ(z) != t_VEC && typ(z) != t_COL) || glength(z) != m)
        pari_err_DIM("my_triple_massey_contract [character]");

    long inputs = lg(gel(T, 1)) - 1;
    GEN result = zerocol(inputs);
    for (long i = 1; i <= m; ++i)
        for (long j = 1; j <= m; ++j)
            for (long k = 1; k <= m; ++k)
            {
                GEN coefficient =
                    Fp_mul(
                        Fp_mul(modii(gel(x, i), p),
                               modii(gel(y, j), p), p),
                        modii(gel(z, k), p), p);
                if (signe(coefficient))
                {
                    GEN word =
                        gel(T, word_index(m, i, j, k));
                    for (long l = 1; l <= inputs; ++l)
                        gel(result, l) =
                            Fp_add(
                                gel(result, l),
                                Fp_mul(coefficient, gel(word, l), p),
                                p);
                }
            }
    return gerepilecopy(av, result);
}

void
my_validate_triple_massey_identities(
    GEN p, GEN T, long m)
{
    if (typ(T) != t_MAT || lg(T) != m * m * m + 1)
        tensor_error("word matrix has the wrong number of columns");

    for (long i = 1; i <= m; ++i)
        for (long j = 1; j <= m; ++j)
            for (long k = 1; k <= m; ++k)
            {
                GEN Mijk = gel(T, word_index(m, i, j, k));
                GEN Mkji = gel(T, word_index(m, k, j, i));
                if (!gequal(Mijk, Mkji))
                    tensor_error("outer-symmetry identity failed");

                GEN cyclic =
                    FpV_add(
                        FpV_add(
                            Mijk,
                            gel(T, word_index(m, j, k, i)), p),
                        gel(T, word_index(m, k, i, j)), p);
                if (!gequal0(cyclic))
                    tensor_error("cyclic shuffle identity failed");

                if (i == j && j == k && !equaliu(p, 3)
                    && !gequal0(Mijk))
                    tensor_error("diagonal triple identity failed");
            }
}
