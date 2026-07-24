// MIT License

/**
 * @file relation_degree3.c
 * @brief Degree-three relation fixtures and Anick strong-freeness tests.
 *
 * The module works entirely over F_5.  It constructs the verified fixture,
 * performs ordered row reduction for candidate monomial orders, and checks
 * combinatorial freeness of the resulting leading words.
 */

#include <string.h>
#include <stdlib.h>
#include <pari/pari.h>
#include "../headers/relation_degree3.h"

#define REL_P 5
#define REL_ROWS 3
#define REL_WORDS 27
#define REL_LIE_DIM 8

static int rel_current_p = REL_P;

static const char rel_letters[] = "abc";
static const int rel_permutations[6][3] = {
    {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
    {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
};

static int
rel_mod(int value)
{
    value %= rel_current_p;
    return value < 0 ? value + rel_current_p : value;
}

static int
rel_inverse(int value)
{
    value = rel_mod(value);
    for (int candidate = 1; candidate < rel_current_p; ++candidate)
        if (rel_mod(value * candidate) == 1) return candidate;
    pari_err(e_MISC, "relation_degree3: attempted to invert zero");
    return 0;
}

static void
rel_decode_word(int index, int word[3])
{
    word[0] = index / 9;
    word[1] = (index / 3) % 3;
    word[2] = index % 3;
}

static int
rel_encode_word(int x, int y, int z)
{
    return 9 * x + 3 * y + z;
}

static void
rel_word_string(int index, char text[4])
{
    int word[3];
    rel_decode_word(index, word);
    for (int i = 0; i < 3; ++i) text[i] = rel_letters[word[i]];
    text[3] = '\0';
}

static int
rel_is_subword(const int *small, int small_length,
               const int *large, int large_length)
{
    if (small_length > large_length) return 0;
    for (int start = 0; start <= large_length - small_length; ++start)
    {
        int equal = 1;
        for (int i = 0; i < small_length; ++i)
            if (small[i] != large[start + i]) { equal = 0; break; }
        if (equal) return 1;
    }
    return 0;
}

static int
rel_words_combinatorially_free(
    const int *words, const int *lengths, int count)
{
    int decoded[REL_WORDS][3];
    for (int i = 0; i < count; ++i)
        rel_decode_word(words[i], decoded[i]);

    for (int i = 0; i < count; ++i)
        for (int j = 0; j < count; ++j)
        {
            if (i != j
                && rel_is_subword(
                    decoded[i], lengths[i], decoded[j], lengths[j]))
                return 0;
            int maximum =
                lengths[i] < lengths[j] ? lengths[i] : lengths[j];
            for (int overlap = 1; overlap < maximum; ++overlap)
            {
                int equal = 1;
                for (int q = 0; q < overlap; ++q)
                    if (decoded[i][q]
                        != decoded[j][lengths[j] - overlap + q])
                    {
                        equal = 0;
                        break;
                    }
                if (equal) return 0;
            }
        }
    return 1;
}

static int
rel_cf_degree3(const int *words, int count)
{
    int lengths[REL_WORDS];
    for (int i = 0; i < count; ++i) lengths[i] = 3;
    return rel_words_combinatorially_free(words, lengths, count);
}

static int
rel_matrix_rank(int rows, int columns, int matrix[REL_ROWS][REL_WORDS])
{
    int copy[REL_ROWS][REL_WORDS];
    memcpy(copy, matrix, sizeof(copy));
    int rank = 0;
    for (int column = 0; column < columns && rank < rows; ++column)
    {
        int pivot = rank;
        while (pivot < rows && !copy[pivot][column]) ++pivot;
        if (pivot == rows) continue;
        if (pivot != rank)
            for (int j = 0; j < columns; ++j)
            {
                int swap = copy[rank][j];
                copy[rank][j] = copy[pivot][j];
                copy[pivot][j] = swap;
            }
        int inverse = rel_inverse(copy[rank][column]);
        for (int j = 0; j < columns; ++j)
            copy[rank][j] = rel_mod(copy[rank][j] * inverse);
        for (int i = 0; i < rows; ++i)
            if (i != rank && copy[i][column])
            {
                int factor = copy[i][column];
                for (int j = 0; j < columns; ++j)
                    copy[i][j] =
                        rel_mod(copy[i][j] - factor * copy[rank][j]);
            }
        ++rank;
    }
    return rank;
}

static int
rel_column_rank(
    int T[REL_ROWS][REL_WORDS], const int *columns, int count)
{
    int selected[REL_ROWS][REL_WORDS] = {{0}};
    for (int row = 0; row < REL_ROWS; ++row)
        for (int j = 0; j < count; ++j)
            selected[row][j] = T[row][columns[j]];
    return rel_matrix_rank(REL_ROWS, count, selected);
}

static int
rel_word_compare(
    int left, int right, const int positions[3], const int weights[3])
{
    int x[3], y[3], x_weight = 0, y_weight = 0;
    rel_decode_word(left, x);
    rel_decode_word(right, y);
    for (int i = 0; i < 3; ++i)
    {
        x_weight += weights[x[i]];
        y_weight += weights[y[i]];
    }
    if (x_weight != y_weight) return x_weight > y_weight ? 1 : -1;
    for (int i = 0; i < 3; ++i)
        if (positions[x[i]] != positions[y[i]])
            return positions[x[i]] > positions[y[i]] ? 1 : -1;
    return 0;
}

static void
rel_make_order(
    const int permutation[3], const int weights[3], int order[REL_WORDS])
{
    int positions[3];
    for (int i = 0; i < 3; ++i) positions[permutation[i]] = i;
    for (int i = 0; i < REL_WORDS; ++i) order[i] = i;
    for (int i = 0; i < REL_WORDS; ++i)
        for (int j = i + 1; j < REL_WORDS; ++j)
            if (rel_word_compare(order[j], order[i], positions, weights) > 0)
            {
                int swap = order[i];
                order[i] = order[j];
                order[j] = swap;
            }
}

static void
rel_ordered_rref(
    int original[REL_ROWS][REL_WORDS], const int order[REL_WORDS],
    int reduced[REL_ROWS][REL_WORDS], int operations[REL_ROWS][REL_ROWS],
    int pivots[REL_ROWS])
{
    memcpy(reduced, original, sizeof(int) * REL_ROWS * REL_WORDS);
    memset(operations, 0, sizeof(int) * REL_ROWS * REL_ROWS);
    for (int i = 0; i < REL_ROWS; ++i) operations[i][i] = 1;

    int rank = 0;
    for (int q = 0; q < REL_WORDS && rank < REL_ROWS; ++q)
    {
        int column = order[q], pivot = rank;
        while (pivot < REL_ROWS && !reduced[pivot][column]) ++pivot;
        if (pivot == REL_ROWS) continue;
        if (pivot != rank)
        {
            for (int j = 0; j < REL_WORDS; ++j)
            {
                int swap = reduced[rank][j];
                reduced[rank][j] = reduced[pivot][j];
                reduced[pivot][j] = swap;
            }
            for (int j = 0; j < REL_ROWS; ++j)
            {
                int swap = operations[rank][j];
                operations[rank][j] = operations[pivot][j];
                operations[pivot][j] = swap;
            }
        }
        int inverse = rel_inverse(reduced[rank][column]);
        for (int j = 0; j < REL_WORDS; ++j)
            reduced[rank][j] = rel_mod(reduced[rank][j] * inverse);
        for (int j = 0; j < REL_ROWS; ++j)
            operations[rank][j] = rel_mod(operations[rank][j] * inverse);
        for (int i = 0; i < REL_ROWS; ++i)
            if (i != rank && reduced[i][column])
            {
                int factor = reduced[i][column];
                for (int j = 0; j < REL_WORDS; ++j)
                    reduced[i][j] =
                        rel_mod(reduced[i][j] - factor * reduced[rank][j]);
                for (int j = 0; j < REL_ROWS; ++j)
                    operations[i][j] =
                        rel_mod(
                            operations[i][j]
                            - factor * operations[rank][j]);
            }
        pivots[rank++] = column;
    }
    if (rank != REL_ROWS)
        pari_err(e_MISC, "relation_degree3: fixture row rank is not 3");

    for (int row = 0; row < REL_ROWS; ++row)
    {
        int found = -1;
        for (int q = 0; q < REL_WORDS; ++q)
            if (reduced[row][order[q]]) { found = order[q]; break; }
        if (found != pivots[row])
            pari_err(e_MISC, "relation_degree3: reported pivot is not leading");
    }
}

static void
rel_print_int_matrix(
    const char *label, int rows, int columns,
    const int *matrix, int stride)
{
    pari_printf("%s[", label);
    for (int i = 0; i < rows; ++i)
    {
        if (i) pari_printf("; ");
        for (int j = 0; j < columns; ++j)
        {
            if (j) pari_printf(",");
            pari_printf("%d", matrix[i * stride + j]);
        }
    }
    pari_printf("]\n");
}

static void
rel_build_lie_basis(int L[REL_LIE_DIM][REL_WORDS])
{
    const int triples[REL_LIE_DIM][3] = {
        {0,1,0}, {0,1,1}, {0,1,2}, {0,2,0},
        {0,2,1}, {0,2,2}, {1,2,1}, {1,2,2}
    };
    memset(L, 0, sizeof(int) * REL_LIE_DIM * REL_WORDS);
    for (int q = 0; q < REL_LIE_DIM; ++q)
    {
        int x = triples[q][0], y = triples[q][1], z = triples[q][2];
        L[q][rel_encode_word(x,y,z)] =
            rel_mod(L[q][rel_encode_word(x,y,z)] + 1);
        L[q][rel_encode_word(y,x,z)] =
            rel_mod(L[q][rel_encode_word(y,x,z)] - 1);
        L[q][rel_encode_word(z,x,y)] =
            rel_mod(L[q][rel_encode_word(z,x,y)] - 1);
        L[q][rel_encode_word(z,y,x)] =
            rel_mod(L[q][rel_encode_word(z,y,x)] + 1);
    }
}

static void
rel_lie_coordinates(
    int T[REL_ROWS][REL_WORDS], int L[REL_LIE_DIM][REL_WORDS],
    int coordinates[REL_ROWS][REL_LIE_DIM])
{
    for (int target = 0; target < REL_ROWS; ++target)
    {
        int augmented[REL_WORDS][REL_LIE_DIM + 1];
        for (int equation = 0; equation < REL_WORDS; ++equation)
        {
            for (int q = 0; q < REL_LIE_DIM; ++q)
                augmented[equation][q] = L[q][equation];
            augmented[equation][REL_LIE_DIM] = T[target][equation];
        }
        int rank = 0, pivot_rows[REL_LIE_DIM];
        for (int column = 0; column < REL_LIE_DIM; ++column)
        {
            int pivot = rank;
            while (pivot < REL_WORDS && !augmented[pivot][column]) ++pivot;
            if (pivot == REL_WORDS)
                pari_err(e_MISC, "relation_degree3: Lie basis rank is not 8");
            for (int j = column; j <= REL_LIE_DIM; ++j)
            {
                int swap = augmented[rank][j];
                augmented[rank][j] = augmented[pivot][j];
                augmented[pivot][j] = swap;
            }
            int inverse = rel_inverse(augmented[rank][column]);
            for (int j = column; j <= REL_LIE_DIM; ++j)
                augmented[rank][j] =
                    rel_mod(augmented[rank][j] * inverse);
            for (int i = 0; i < REL_WORDS; ++i)
                if (i != rank && augmented[i][column])
                {
                    int factor = augmented[i][column];
                    for (int j = column; j <= REL_LIE_DIM; ++j)
                        augmented[i][j] =
                            rel_mod(
                                augmented[i][j]
                                - factor * augmented[rank][j]);
                }
            pivot_rows[column] = rank++;
        }
        for (int i = rank; i < REL_WORDS; ++i)
            if (augmented[i][REL_LIE_DIM])
                pari_err(e_MISC, "relation_degree3: row is outside free Lie L3");
        for (int q = 0; q < REL_LIE_DIM; ++q)
            coordinates[target][q] =
                augmented[pivot_rows[q]][REL_LIE_DIM];

        for (int column = 0; column < REL_WORDS; ++column)
        {
            int reconstructed = 0;
            for (int q = 0; q < REL_LIE_DIM; ++q)
                reconstructed += coordinates[target][q] * L[q][column];
            if (rel_mod(reconstructed) != T[target][column])
                pari_err(e_MISC, "relation_degree3: Lie reconstruction failed");
        }
    }
}

static int
rel_test_anick(
    int T[REL_ROWS][REL_WORDS], const int order[REL_WORDS],
    int verbose, int print_certificate)
{
    int reduced[REL_ROWS][REL_WORDS], operations[REL_ROWS][REL_ROWS];
    int pivots[REL_ROWS];
    rel_ordered_rref(T, order, reduced, operations, pivots);
    int free = rel_cf_degree3(pivots, REL_ROWS);
    char words[REL_ROWS][4];
    for (int i = 0; i < REL_ROWS; ++i)
        rel_word_string(pivots[i], words[i]);
    if (verbose)
        pari_printf(
            "    Anick pivots = [%s,%s,%s], combinatorially free = %s\n",
            words[0], words[1], words[2], free ? "yes" : "no");
    if (free || print_certificate)
    {
        rel_print_int_matrix(
            "    row-operation matrix U = ", REL_ROWS, REL_ROWS,
            &operations[0][0], REL_ROWS);
        rel_print_int_matrix(
            "    transformed relation matrix U*T = ",
            REL_ROWS, REL_WORDS, &reduced[0][0], REL_WORDS);
    }
    return free;
}

static int
rel_test_efrat(
    int T[REL_ROWS][REL_WORDS], const int order[REL_WORDS],
    int verbose, int print_candidates)
{
    int support[REL_WORDS], support_count = 0;
    for (int q = 0; q < REL_WORDS; ++q)
    {
        int column = order[q];
        if (T[0][column] || T[1][column] || T[2][column])
            support[support_count++] = column;
    }
    if (verbose)
        pari_printf("    Efrat nonzero support size = %d\n", support_count);
    int successes = 0;
    for (int count = 1; count <= support_count; ++count)
    {
        int rank = rel_column_rank(T, support, count);
        if (rank != REL_ROWS) continue;
        int free = rel_cf_degree3(support, count);
        if (print_candidates)
            pari_printf(
                "      upper tail size=%d rank=3 CF=%s\n",
                count, free ? "yes" : "no");
        if (free)
        {
            ++successes;
            pari_printf("      Efrat certificate B = [");
            for (int i = 0; i < count; ++i)
            {
                char word[4];
                rel_word_string(support[i], word);
                pari_printf("%s%s", i ? "," : "", word);
            }
            pari_printf("]\n");
            int alpha[REL_ROWS][REL_WORDS] = {{0}};
            for (int row = 0; row < REL_ROWS; ++row)
                for (int i = 0; i < count; ++i)
                    alpha[row][i] = T[row][support[i]];
            rel_print_int_matrix(
                "      alpha_B = ", REL_ROWS, count,
                &alpha[0][0], REL_WORDS);
            pari_printf("      rank(alpha_B) = 3\n");
        }
    }
    if (verbose)
        pari_printf(
            "    Efrat successful upper tails = %d\n", successes);
    return successes != 0;
}

static void
rel_combinatorial_unit_tests(void)
{
    int pass[] = {
        rel_encode_word(0,0,1),
        rel_encode_word(0,2,1),
        rel_encode_word(2,0,1)
    };
    int overlap_fail[] = {
        rel_encode_word(0,1,2),
        rel_encode_word(1,2,0)
    };
    int duplicate_fail[] = {
        rel_encode_word(0,1,2),
        rel_encode_word(0,1,2)
    };
    if (!rel_cf_degree3(pass, 3)
        || rel_cf_degree3(overlap_fail, 2)
        || rel_cf_degree3(duplicate_fail, 2))
        pari_err(e_MISC, "relation_degree3: combinatorial-freeness unit test failed");
    pari_printf(
        "  CF unit tests: [aab,acb,cab] passes; "
        "[abc,bca] overlap fails; duplicate [abc,abc] fails\n");
}

static void
rel_census(int T[REL_ROWS][REL_WORDS])
{
    int count = 0, printed = 0;
    pari_printf("  independent combinatorially-free triple sample: ");
    for (int i = 0; i < REL_WORDS; ++i)
        for (int j = i + 1; j < REL_WORDS; ++j)
            for (int k = j + 1; k < REL_WORDS; ++k)
            {
                int triple[3] = {i, j, k};
                if (!rel_cf_degree3(triple, 3)
                    || rel_column_rank(T, triple, 3) != 3)
                    continue;
                ++count;
                if (printed < 10)
                {
                    char wi[4], wj[4], wk[4];
                    rel_word_string(i, wi);
                    rel_word_string(j, wj);
                    rel_word_string(k, wk);
                    pari_printf(
                        "%s[%s,%s,%s]", printed ? ", " : "",
                        wi, wj, wk);
                    ++printed;
                }
            }
    pari_printf("\n  census total = %d\n", count);
    pari_printf(
        "  THIS ALONE IS NOT A STRONG-FREENESS CERTIFICATE.\n");
}

static int
rel_gcd(int a, int b)
{
    while (b) { int remainder = a % b; a = b; b = remainder; }
    return a;
}

static void
rel_weighted_search(int T[REL_ROWS][REL_WORDS])
{
    int seen[500][REL_WORDS], seen_count = 0;
    int any = 0;
    for (int wa = 1; wa <= 8; ++wa)
        for (int wb = 1; wb <= 8; ++wb)
            for (int wc = 1; wc <= 8; ++wc)
            {
                if (rel_gcd(rel_gcd(wa, wb), wc) != 1) continue;
                int weights[3] = {wa, wb, wc};
                for (int permutation = 0; permutation < 6; ++permutation)
                {
                    int order[REL_WORDS];
                    rel_make_order(
                        rel_permutations[permutation], weights, order);
                    int duplicate = 0;
                    for (int q = 0; q < seen_count; ++q)
                        if (!memcmp(seen[q], order, sizeof(order)))
                        {
                            duplicate = 1;
                            break;
                        }
                    if (duplicate) continue;
                    memcpy(seen[seen_count++], order, sizeof(order));
                    int anick = rel_test_anick(T, order, 0, 0);
                    int efrat = rel_test_efrat(T, order, 0, 0);
                    if (anick || efrat)
                    {
                        pari_printf(
                            "  weighted certificate: weights=(%d,%d,%d), "
                            "tie order=%c<%c<%c\n",
                            wa, wb, wc,
                            rel_letters[rel_permutations[permutation][0]],
                            rel_letters[rel_permutations[permutation][1]],
                            rel_letters[rel_permutations[permutation][2]]);
                        any = 1;
                    }
                }
            }
    pari_printf(
        "  distinct weighted degree-3 orders tested = %d\n", seen_count);
    if (!any)
        pari_printf(
            "  no certificate was found in the tested weighted-order family.\n");
}

typedef struct {
    int length;
    int letter[5];
} RelFiniteWord;

typedef struct {
    int value[3];
} RelSigma;

typedef struct {
    int sigma_count;
    RelSigma sigma[3];
    int permutation[3];
} RelSection8Order;

typedef struct {
    int words[3];
    int top_three_count;
    int maximum_pivot_intersection;
    int minimum_extra_above;
} RelPromisingTriple;

static long
rel_sigma_star(const RelFiniteWord *word, const RelSigma *sigma)
{
    long result = 0;
    for (int i = 0; i < word->length; ++i)
        result += sigma->value[word->letter[i]];
    return result;
}

static long
rel_sigma_sharp(const RelFiniteWord *word, const RelSigma *sigma)
{
    long result = 0;
    for (int i = 0; i < word->length; ++i)
        result += (long)(i + 1) * sigma->value[word->letter[i]];
    return result;
}

static RelFiniteWord
rel_finite_word(int length, int code)
{
    RelFiniteWord word;
    word.length = length;
    memset(word.letter, 0, sizeof(word.letter));
    for (int i = length - 1; i >= 0; --i)
    {
        word.letter[i] = code % 3;
        code /= 3;
    }
    return word;
}

static RelFiniteWord
rel_concatenate(const RelFiniteWord *left, const RelFiniteWord *right)
{
    RelFiniteWord result;
    result.length = left->length + right->length;
    if (result.length > 5)
        pari_err(e_MISC, "relation_degree3: finite test word is too long");
    for (int i = 0; i < left->length; ++i)
        result.letter[i] = left->letter[i];
    for (int i = 0; i < right->length; ++i)
        result.letter[left->length + i] = right->letter[i];
    return result;
}

static int
rel_section8_compare(
    const RelFiniteWord *left, const RelFiniteWord *right,
    const RelSection8Order *order)
{
    if (left->length != right->length)
        return left->length < right->length ? -1 : 1;
    for (int q = order->sigma_count - 1; q >= 0; --q)
    {
        long left_star = rel_sigma_star(left, &order->sigma[q]);
        long right_star = rel_sigma_star(right, &order->sigma[q]);
        if (left_star != right_star)
            return left_star < right_star ? -1 : 1;
        long left_sharp = rel_sigma_sharp(left, &order->sigma[q]);
        long right_sharp = rel_sigma_sharp(right, &order->sigma[q]);
        if (left_sharp != right_sharp)
            return left_sharp < right_sharp ? -1 : 1;
    }
    int positions[3];
    for (int i = 0; i < 3; ++i)
        positions[order->permutation[i]] = i;
    for (int i = 0; i < left->length; ++i)
        if (positions[left->letter[i]] != positions[right->letter[i]])
            return positions[left->letter[i]] < positions[right->letter[i]]
                ? -1 : 1;
    return 0;
}

static void
rel_make_section8_order(
    const RelSection8Order *section8, int ordered_words[REL_WORDS])
{
    RelFiniteWord words[REL_WORDS];
    for (int i = 0; i < REL_WORDS; ++i)
    {
        ordered_words[i] = i;
        words[i] = rel_finite_word(3, i);
    }
    for (int i = 0; i < REL_WORDS; ++i)
        for (int j = i + 1; j < REL_WORDS; ++j)
            if (rel_section8_compare(
                    &words[ordered_words[j]],
                    &words[ordered_words[i]], section8) > 0)
            {
                int swap = ordered_words[i];
                ordered_words[i] = ordered_words[j];
                ordered_words[j] = swap;
            }
}

static void
rel_section8_unit_tests(void)
{
    RelSigma sigma = {{1, 2, 4}};
    RelFiniteWord abc = {0}, cba = {0};
    abc.length = cba.length = 3;
    abc.letter[0] = 0; abc.letter[1] = 1; abc.letter[2] = 2;
    cba.letter[0] = 2; cba.letter[1] = 1; cba.letter[2] = 0;
    if (rel_sigma_star(&abc, &sigma) != 7
        || rel_sigma_star(&cba, &sigma) != 7
        || rel_sigma_sharp(&abc, &sigma) != 17
        || rel_sigma_sharp(&cba, &sigma) != 11)
        pari_err(e_MISC, "relation_degree3: sigma star/sharp unit test failed");

    const RelSigma tests[] = {
        {{1,2,4}}, {{1,0,1}}, {{0,3,2}}
    };
    for (int test = 0; test < 3; ++test)
        for (int left_length = 0; left_length <= 4; ++left_length)
        {
            int left_count = 1;
            for (int i = 0; i < left_length; ++i) left_count *= 3;
            for (int right_length = 0;
                 right_length + left_length <= 4; ++right_length)
            {
                int right_count = 1;
                for (int i = 0; i < right_length; ++i) right_count *= 3;
                for (int left_code = 0; left_code < left_count; ++left_code)
                    for (int right_code = 0;
                         right_code < right_count; ++right_code)
                    {
                        RelFiniteWord left =
                            rel_finite_word(left_length, left_code);
                        RelFiniteWord right =
                            rel_finite_word(right_length, right_code);
                        RelFiniteWord joined =
                            rel_concatenate(&left, &right);
                        long expected_star =
                            rel_sigma_star(&left, &tests[test])
                            + rel_sigma_star(&right, &tests[test]);
                        long expected_sharp =
                            (long)left.length
                                * rel_sigma_star(&right, &tests[test])
                            + rel_sigma_sharp(&left, &tests[test])
                            + rel_sigma_sharp(&right, &tests[test]);
                        if (rel_sigma_star(&joined, &tests[test])
                                != expected_star
                            || rel_sigma_sharp(&joined, &tests[test])
                                != expected_sharp)
                            pari_err(
                                e_MISC,
                                "relation_degree3: sigma concatenation "
                                "identity failed");
                    }
            }
        }
    pari_printf(
        "\n  Efrat 8.1 units: sigma=[1,2,4], "
        "star(abc)=star(cba)=7, sharp(abc)=17, sharp(cba)=11\n");
    pari_printf(
        "  exhaustive star/sharp concatenation identities through length 4 = OK\n");
}

static void
rel_section8_monoid_checks(void)
{
    const RelSigma representative[] = {
        {{1,2,4}}, {{1,0,1}}
    };
    int comparators = 0;
    for (int sigma_choice = 0; sigma_choice < 2; ++sigma_choice)
        for (int permutation = 0; permutation < 6; ++permutation)
        {
            RelSection8Order order;
            order.sigma_count = 1;
            order.sigma[0] = representative[sigma_choice];
            memcpy(
                order.permutation, rel_permutations[permutation],
                sizeof(order.permutation));
            RelFiniteWord empty = rel_finite_word(0, 0);
            for (int length = 1; length <= 3; ++length)
            {
                int count = 1;
                for (int i = 0; i < length; ++i) count *= 3;
                for (int code = 0; code < count; ++code)
                {
                    RelFiniteWord word = rel_finite_word(length, code);
                    if (rel_section8_compare(&empty, &word, &order) >= 0)
                        pari_err(e_MISC, "relation_degree3: empty word order failed");
                }
            }
            for (int left_length = 0; left_length <= 3; ++left_length)
            {
                int left_count = 1;
                for (int i = 0; i < left_length; ++i) left_count *= 3;
                for (int right_length = 0; right_length <= 3; ++right_length)
                {
                    int right_count = 1;
                    for (int i = 0; i < right_length; ++i) right_count *= 3;
                    for (int left_code = 0; left_code < left_count; ++left_code)
                        for (int right_code = 0;
                             right_code < right_count; ++right_code)
                        {
                            RelFiniteWord left =
                                rel_finite_word(left_length, left_code);
                            RelFiniteWord right =
                                rel_finite_word(right_length, right_code);
                            if (rel_section8_compare(&left, &right, &order) >= 0)
                                continue;
                            for (int u_length = 0; u_length <= 1; ++u_length)
                                for (int v_length = 0; v_length <= 1; ++v_length)
                                {
                                    if (u_length + left_length + v_length > 5
                                        || u_length + right_length + v_length > 5)
                                        continue;
                                    int u_count = u_length ? 3 : 1;
                                    int v_count = v_length ? 3 : 1;
                                    for (int u_code = 0; u_code < u_count; ++u_code)
                                        for (int v_code = 0;
                                             v_code < v_count; ++v_code)
                                        {
                                            RelFiniteWord u =
                                                rel_finite_word(u_length, u_code);
                                            RelFiniteWord v =
                                                rel_finite_word(v_length, v_code);
                                            RelFiniteWord ul =
                                                rel_concatenate(&u, &left);
                                            RelFiniteWord ur =
                                                rel_concatenate(&u, &right);
                                            RelFiniteWord ulv =
                                                rel_concatenate(&ul, &v);
                                            RelFiniteWord urv =
                                                rel_concatenate(&ur, &v);
                                            if (rel_section8_compare(
                                                    &ulv, &urv, &order) >= 0)
                                                pari_err(
                                                    e_MISC,
                                                    "relation_degree3: ordered-"
                                                    "monoid compatibility failed");
                                        }
                                }
                        }
                }
            }
            ++comparators;
        }
    pari_printf(
        "  ordered-monoid implementation checks for %d representative "
        "comparators and words through length 5 = OK\n", comparators);
}

static int
rel_collect_promising(
    int T[REL_ROWS][REL_WORDS], RelPromisingTriple triples[30])
{
    int count = 0;
    for (int i = 0; i < REL_WORDS; ++i)
        for (int j = i + 1; j < REL_WORDS; ++j)
            for (int k = j + 1; k < REL_WORDS; ++k)
            {
                int candidate[3] = {i, j, k};
                if (!rel_cf_degree3(candidate, 3)
                    || rel_column_rank(T, candidate, 3) != 3)
                    continue;
                memcpy(triples[count].words, candidate, sizeof(candidate));
                triples[count].top_three_count = 0;
                triples[count].maximum_pivot_intersection = 0;
                triples[count].minimum_extra_above = REL_WORDS;
                ++count;
            }
    return count;
}

static int
rel_same_triple(const int left[3], const int right[3])
{
    for (int i = 0; i < 3; ++i)
    {
        int found = 0;
        for (int j = 0; j < 3; ++j)
            if (left[i] == right[j]) { found = 1; break; }
        if (!found) return 0;
    }
    return 1;
}

static void
rel_track_promising(
    int T[REL_ROWS][REL_WORDS], const int order[REL_WORDS],
    RelPromisingTriple *triples, int triple_count)
{
    int reduced[REL_ROWS][REL_WORDS], operations[REL_ROWS][REL_ROWS];
    int pivots[3], support[REL_WORDS], support_count = 0;
    rel_ordered_rref(T, order, reduced, operations, pivots);
    for (int q = 0; q < REL_WORDS; ++q)
    {
        int column = order[q];
        if (T[0][column] || T[1][column] || T[2][column])
            support[support_count++] = column;
    }
    for (int t = 0; t < triple_count; ++t)
    {
        if (rel_same_triple(triples[t].words, support))
            ++triples[t].top_three_count;
        int intersection = 0, maximum_position = -1;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
                if (triples[t].words[i] == pivots[j]) ++intersection;
            for (int position = 0; position < support_count; ++position)
                if (triples[t].words[i] == support[position]
                    && position > maximum_position)
                    maximum_position = position;
        }
        if (intersection > triples[t].maximum_pivot_intersection)
            triples[t].maximum_pivot_intersection = intersection;
        int extra = maximum_position + 1 - 3;
        if (extra < triples[t].minimum_extra_above)
            triples[t].minimum_extra_above = extra;
    }
}

static int
rel_is_partition_sequence(const RelSection8Order *order)
{
    int used[3] = {0, 0, 0};
    for (int q = 0; q < order->sigma_count; ++q)
    {
        int nonempty = 0;
        for (int letter = 0; letter < 3; ++letter)
        {
            int value = order->sigma[q].value[letter];
            if (value != 0 && value != 1) return 0;
            if (value)
            {
                if (used[letter]) return 0;
                used[letter] = 1;
                nonempty = 1;
            }
        }
        if (!nonempty) return 0;
    }
    int y0_nonempty = 0;
    for (int letter = 0; letter < 3; ++letter)
        if (!used[letter]) y0_nonempty = 1;
    return y0_nonempty;
}

static int
rel_section8_efrat_success(
    int T[REL_ROWS][REL_WORDS], const int order[REL_WORDS])
{
    int support[REL_WORDS], count = 0;
    for (int q = 0; q < REL_WORDS; ++q)
    {
        int column = order[q];
        if (T[0][column] || T[1][column] || T[2][column])
            support[count++] = column;
    }
    for (int size = 1; size <= count; ++size)
        if (rel_column_rank(T, support, size) == 3
            && rel_cf_degree3(support, size))
            return 1;
    return 0;
}

static long
rel_integer_power(long base, int exponent)
{
    long result = 1;
    while (exponent--) result *= base;
    return result;
}

static void
rel_search_section8_stage(
    int T[REL_ROWS][REL_WORDS], const RelSigma *maps, int map_count,
    int sequence_length, const char *stage,
    RelPromisingTriple *triples, int triple_count,
    int *total_anick, int *total_efrat)
{
    int seen[200][REL_WORDS], seen_count = 0;
    int seen_anick[200] = {0}, seen_efrat[200] = {0};
    int seen_partition[200] = {0};
    int anick_count = 0, efrat_count = 0;
    long sequence_count = rel_integer_power(map_count, sequence_length);
    for (long code = 0; code < sequence_count; ++code)
    {
        long remainder = code;
        RelSection8Order base;
        base.sigma_count = sequence_length;
        for (int q = 0; q < sequence_length; ++q)
        {
            base.sigma[q] = maps[remainder % map_count];
            remainder /= map_count;
        }
        for (int permutation = 0; permutation < 6; ++permutation)
        {
            memcpy(
                base.permutation, rel_permutations[permutation],
                sizeof(base.permutation));
            int order[REL_WORDS];
            rel_make_section8_order(&base, order);
            int duplicate = -1;
            for (int q = 0; q < seen_count; ++q)
                if (!memcmp(seen[q], order, sizeof(order)))
                {
                    duplicate = q;
                    break;
                }
            if (duplicate >= 0)
            {
                if (rel_is_partition_sequence(&base))
                    seen_partition[duplicate] = 1;
                continue;
            }
            int order_number = seen_count++;
            memcpy(seen[order_number], order, sizeof(order));
            int reduced[REL_ROWS][REL_WORDS], operations[REL_ROWS][REL_ROWS];
            int pivots[3];
            rel_ordered_rref(T, order, reduced, operations, pivots);
            int anick = rel_cf_degree3(pivots, 3);
            int efrat = rel_section8_efrat_success(T, order);
            seen_anick[order_number] = anick;
            seen_efrat[order_number] = efrat;
            seen_partition[order_number] =
                rel_is_partition_sequence(&base);
            anick_count += anick;
            efrat_count += efrat;
            rel_track_promising(T, order, triples, triple_count);
            if (anick || efrat)
            {
                pari_printf(
                    "  CERTIFICATE order: %s s=%d sigmas=",
                    stage, sequence_length);
                for (int q = 0; q < sequence_length; ++q)
                    pari_printf(
                        "%s[%d,%d,%d]", q ? "," : "",
                        base.sigma[q].value[0],
                        base.sigma[q].value[1],
                        base.sigma[q].value[2]);
                pari_printf(
                    " tie=%c<%c<%c Anick=%s Efrat=%s\n",
                    rel_letters[base.permutation[0]],
                    rel_letters[base.permutation[1]],
                    rel_letters[base.permutation[2]],
                    anick ? "yes" : "no", efrat ? "yes" : "no");
                rel_test_anick(T, order, 1, anick);
                rel_test_efrat(T, order, 1, efrat);
            }
        }
    }
    long raw = sequence_count * 6;
    int partition_orders = 0, partition_anick = 0, partition_efrat = 0;
    for (int q = 0; q < seen_count; ++q)
        if (seen_partition[q])
        {
            ++partition_orders;
            partition_anick += seen_anick[q];
            partition_efrat += seen_efrat[q];
        }
    pari_printf(
        "  %s s=%d: raw=%ld distinct=%d Anick=%d Efrat=%d",
        stage, sequence_length, raw, seen_count,
        anick_count, efrat_count);
    if (!strcmp(stage, "Stage1"))
        pari_printf(
            " partition-derived-orders=%d partition-Anick=%d "
            "partition-Efrat=%d",
            partition_orders, partition_anick, partition_efrat);
    pari_printf("\n");
    *total_anick += anick_count;
    *total_efrat += efrat_count;
}

static int
rel_normalized_sigma_maps(RelSigma maps[125])
{
    int count = 0;
    for (int a = 0; a <= 4; ++a)
        for (int b = 0; b <= 4; ++b)
            for (int c = 0; c <= 4; ++c)
            {
                int minimum = a;
                if (b < minimum) minimum = b;
                if (c < minimum) minimum = c;
                int values[3] = {a - minimum, b - minimum, c - minimum};
                if (!values[0] && !values[1] && !values[2]) continue;
                int gcd = rel_gcd(rel_gcd(values[0], values[1]), values[2]);
                for (int i = 0; i < 3; ++i) values[i] /= gcd;
                int duplicate = 0;
                for (int q = 0; q < count; ++q)
                    if (maps[q].value[0] == values[0]
                        && maps[q].value[1] == values[1]
                        && maps[q].value[2] == values[2])
                    {
                        duplicate = 1;
                        break;
                    }
                if (!duplicate)
                {
                    for (int i = 0; i < 3; ++i)
                        maps[count].value[i] = values[i];
                    ++count;
                }
            }
    return count;
}

static void
rel_run_section8_search(int T[REL_ROWS][REL_WORDS])
{
    rel_section8_unit_tests();
    rel_section8_monoid_checks();

    RelPromisingTriple triples[30];
    int triple_count = rel_collect_promising(T, triples);
    if (triple_count != 24)
        pari_err(e_MISC, "relation_degree3: promising triple count changed");

    RelSigma characteristic[6];
    int map_count = 0;
    pari_printf("  Stage1 characteristic maps: ");
    for (int mask = 1; mask < 7; ++mask)
    {
        for (int letter = 0; letter < 3; ++letter)
            characteristic[map_count].value[letter] =
                (mask >> letter) & 1;
        pari_printf(
            "%s[%d,%d,%d]", map_count ? "," : "",
            characteristic[map_count].value[0],
            characteristic[map_count].value[1],
            characteristic[map_count].value[2]);
        ++map_count;
    }
    pari_printf("\n");
    pari_printf(
        "  exact partition sequences: s=1 has 6 ordered (Y0,Y1) "
        "partitions; s=2 has 6 ordered (Y0,Y1,Y2) singleton partitions; "
        "s=3 has 0 with all blocks nonempty\n");

    int stage1_anick = 0, stage1_efrat = 0;
    for (int length = 1; length <= 3; ++length)
        rel_search_section8_stage(
            T, characteristic, map_count, length, "Stage1",
            triples, triple_count, &stage1_anick, &stage1_efrat);

    int stage2_anick = 0, stage2_efrat = 0;
    if (!stage1_anick && !stage1_efrat)
    {
        RelSigma normalized[125];
        int normalized_count = rel_normalized_sigma_maps(normalized);
        pari_printf(
            "  Stage2 normalized sigma maps = %d\n", normalized_count);
        for (int length = 1; length <= 2; ++length)
            rel_search_section8_stage(
                T, normalized, normalized_count, length, "Stage2",
                triples, triple_count, &stage2_anick, &stage2_efrat);
    }

    pari_printf("  behavior of 24 promising triples:\n");
    for (int q = 0; q < triple_count; ++q)
    {
        char first[4], second[4], third[4];
        rel_word_string(triples[q].words[0], first);
        rel_word_string(triples[q].words[1], second);
        rel_word_string(triples[q].words[2], third);
        pari_printf(
            "    [%s,%s,%s]: top3=%d max-pivot-intersection=%d "
            "min-extra-above=%d\n",
            first, second, third, triples[q].top_three_count,
            triples[q].maximum_pivot_intersection,
            triples[q].minimum_extra_above);
    }
    int all_anick = stage1_anick + stage2_anick;
    int all_efrat = stage1_efrat + stage2_efrat;
    pari_printf(
        "  Section8 totals: Anick certificates=%d Efrat certificates=%d\n",
        all_anick, all_efrat);
    if (!all_anick && !all_efrat)
        pari_printf(
            "  No certificate was found in the tested Efrat 8.1 "
            "order families.\n");
}

void
my_run_relation_degree3_fixture(GEN fixture, GEN p)
{
    if (!equaliu(p, REL_P) || typ(fixture) != t_MAT
        || glength(fixture) != REL_WORDS
        || lg(gel(fixture, 1)) != REL_ROWS + 1)
        pari_err(e_MISC, "relation_degree3: fixture must be 3 x 27 over F_5");
    rel_current_p = REL_P;

    int T[REL_ROWS][REL_WORDS];
    for (int column = 0; column < REL_WORDS; ++column)
        for (int row = 0; row < REL_ROWS; ++row)
            T[row][column] =
                (int)umodiu(gcoeff(fixture, row + 1, column + 1), REL_P);
    if (rel_matrix_rank(REL_ROWS, REL_WORDS, T) != REL_ROWS)
        pari_err(e_MISC, "relation_degree3: fixture rank is not 3");

    pari_printf("MASSEY_RELATION_TEST pure algebra fixture\n");
    pari_printf("  shape(T) = 3 x 27; rank(T) = 3\n");
    pari_printf(
        "  convention: epsilon_ijk(r_l) = T[l,ijk] in degree 3; "
        "R3 = rowspace(T)\n");
    pari_printf(
        "  generator rank=3; relation rank=3; dim(R3)=3; "
        "Zassenhaus initial degree=3\n");

    int L[REL_LIE_DIM][REL_WORDS];
    rel_build_lie_basis(L);
    /* Check the eight rows with a local 8-row elimination. */
    int lie_copy[REL_LIE_DIM][REL_WORDS];
    memcpy(lie_copy, L, sizeof(L));
    int lie_rank = 0;
    for (int column = 0; column < REL_WORDS; ++column)
    {
        int pivot = lie_rank;
        while (pivot < REL_LIE_DIM && !lie_copy[pivot][column]) ++pivot;
        if (pivot == REL_LIE_DIM) continue;
        for (int j = 0; j < REL_WORDS; ++j)
        {
            int swap = lie_copy[lie_rank][j];
            lie_copy[lie_rank][j] = lie_copy[pivot][j];
            lie_copy[pivot][j] = swap;
        }
        int inverse = rel_inverse(lie_copy[lie_rank][column]);
        for (int j = 0; j < REL_WORDS; ++j)
            lie_copy[lie_rank][j] =
                rel_mod(lie_copy[lie_rank][j] * inverse);
        for (int i = 0; i < REL_LIE_DIM; ++i)
            if (i != lie_rank && lie_copy[i][column])
            {
                int factor = lie_copy[i][column];
                for (int j = 0; j < REL_WORDS; ++j)
                    lie_copy[i][j] =
                        rel_mod(
                            lie_copy[i][j]
                            - factor * lie_copy[lie_rank][j]);
            }
        ++lie_rank;
    }
    if (lie_rank != REL_LIE_DIM)
        pari_err(e_MISC, "relation_degree3: proposed Lie basis is dependent");
    int coordinates[REL_ROWS][REL_LIE_DIM];
    rel_lie_coordinates(T, L, coordinates);
    pari_printf(
        "  Lie basis rank = 8 for "
        "[[a,b],a],[[a,b],b],[[a,b],c],[[a,c],a],"
        "[[a,c],b],[[a,c],c],[[b,c],b],[[b,c],c]\n");
    rel_print_int_matrix(
        "  relation Lie coordinates = ", REL_ROWS, REL_LIE_DIM,
        &coordinates[0][0], REL_LIE_DIM);
    pari_printf("  all three Lie reconstructions = exact\n");

    rel_combinatorial_unit_tests();
    int unit_weights[3] = {1, 1, 1};
    int any_lex_anick = 0, any_lex_efrat = 0;
    pari_printf("\n  six lexicographic orders:\n");
    for (int permutation = 0; permutation < 6; ++permutation)
    {
        int order[REL_WORDS];
        rel_make_order(
            rel_permutations[permutation], unit_weights, order);
        pari_printf(
            "  order %c<%c<%c\n",
            rel_letters[rel_permutations[permutation][0]],
            rel_letters[rel_permutations[permutation][1]],
            rel_letters[rel_permutations[permutation][2]]);
        any_lex_anick |= rel_test_anick(T, order, 1, 0);
        any_lex_efrat |= rel_test_efrat(T, order, 1, 1);
    }

    pari_printf("\n  triple census:\n");
    rel_census(T);

    if (!any_lex_anick && !any_lex_efrat)
    {
        pari_printf("\n  bounded weighted-degree-lex search:\n");
        rel_weighted_search(T);
    }

    if (any_lex_anick)
        pari_printf(
            "  Anick certificate found: initial forms are strongly free.\n");
    if (any_lex_efrat)
        pari_printf(
            "  Efrat Theorem 7.1 certificate found: presentation is mild.\n");
    if (!any_lex_anick && !any_lex_efrat)
        pari_printf(
            "  The tested sufficient criteria did not prove mildness.\n");
    rel_run_section8_search(T);
    pari_printf("MASSEY_RELATION_TEST completed without number-field arithmetic\n");
}

static int
rel_det3(const int matrix[3][3])
{
    return rel_mod(
        matrix[0][0]
            * (matrix[1][1] * matrix[2][2]
               - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
            * (matrix[1][0] * matrix[2][2]
               - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
            * (matrix[1][0] * matrix[2][1]
               - matrix[1][1] * matrix[2][0]));
}

static void
rel_inverse_square(int size, const int *matrix, int *inverse)
{
    int augmented[REL_WORDS][2 * REL_WORDS] = {{0}};
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
        {
            augmented[i][j] = matrix[i * size + j];
            augmented[i][size + j] = i == j;
        }
    for (int column = 0; column < size; ++column)
    {
        int pivot = column;
        while (pivot < size && !augmented[pivot][column]) ++pivot;
        if (pivot == size)
            pari_err(e_MISC, "mild certificate: matrix is singular");
        for (int j = 0; j < 2 * size; ++j)
        {
            int swap = augmented[column][j];
            augmented[column][j] = augmented[pivot][j];
            augmented[pivot][j] = swap;
        }
        int scalar = rel_inverse(augmented[column][column]);
        for (int j = 0; j < 2 * size; ++j)
            augmented[column][j] =
                rel_mod(augmented[column][j] * scalar);
        for (int i = 0; i < size; ++i)
            if (i != column && augmented[i][column])
            {
                int factor = augmented[i][column];
                for (int j = 0; j < 2 * size; ++j)
                    augmented[i][j] =
                        rel_mod(
                            augmented[i][j]
                            - factor * augmented[column][j]);
            }
    }
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            inverse[i * size + j] = augmented[i][size + j];
}

static void
rel_verify_inverse(
    int size, const int *matrix, const int *inverse, const char *name)
{
    for (int side = 0; side < 2; ++side)
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
            {
                int value = 0;
                for (int k = 0; k < size; ++k)
                    value += side
                        ? inverse[i * size + k] * matrix[k * size + j]
                        : matrix[i * size + k] * inverse[k * size + j];
                if (rel_mod(value) != (i == j))
                    pari_err(e_MISC, "mild certificate: inverse check failed");
            }
    pari_printf("  %s inverse identities = OK\n", name);
}

static void
rel_build_cubic_substitution(
    const int degree_one[3][3], int substitution[REL_WORDS][REL_WORDS])
{
    memset(
        substitution, 0,
        sizeof(int) * REL_WORDS * REL_WORDS);
    for (int source = 0; source < REL_WORDS; ++source)
    {
        int word[3];
        rel_decode_word(source, word);
        for (int x = 0; x < 3; ++x)
            for (int y = 0; y < 3; ++y)
                for (int z = 0; z < 3; ++z)
                    substitution[source][rel_encode_word(x, y, z)] =
                        rel_mod(
                            degree_one[x][word[0]]
                            * degree_one[y][word[1]]
                            * degree_one[z][word[2]]);
    }
}

static int
rel_rank27(int matrix[REL_WORDS][REL_WORDS])
{
    int copy[REL_WORDS][REL_WORDS];
    memcpy(copy, matrix, sizeof(copy));
    int rank = 0;
    for (int column = 0; column < REL_WORDS; ++column)
    {
        int pivot = rank;
        while (pivot < REL_WORDS && !copy[pivot][column]) ++pivot;
        if (pivot == REL_WORDS) continue;
        for (int j = 0; j < REL_WORDS; ++j)
        {
            int swap = copy[rank][j];
            copy[rank][j] = copy[pivot][j];
            copy[pivot][j] = swap;
        }
        int scalar = rel_inverse(copy[rank][column]);
        for (int j = 0; j < REL_WORDS; ++j)
            copy[rank][j] = rel_mod(copy[rank][j] * scalar);
        for (int i = 0; i < REL_WORDS; ++i)
            if (i != rank && copy[i][column])
            {
                int factor = copy[i][column];
                for (int j = 0; j < REL_WORDS; ++j)
                    copy[i][j] =
                        rel_mod(copy[i][j] - factor * copy[rank][j]);
            }
        ++rank;
    }
    return rank;
}

static void
rel_apply_cubic_substitution(
    int input[REL_ROWS][REL_WORDS],
    int substitution[REL_WORDS][REL_WORDS],
    int output[REL_ROWS][REL_WORDS])
{
    memset(output, 0, sizeof(int) * REL_ROWS * REL_WORDS);
    for (int row = 0; row < REL_ROWS; ++row)
        for (int source = 0; source < REL_WORDS; ++source)
            for (int target = 0; target < REL_WORDS; ++target)
                output[row][target] =
                    rel_mod(
                        output[row][target]
                        + input[row][source]
                            * substitution[source][target]);
}

static void
rel_compare_matrix(
    const int *actual, const int *expected,
    int rows, int columns, const char *name)
{
    for (int i = 0; i < rows * columns; ++i)
        if (actual[i] != expected[i])
        {
            pari_printf(
                "  %s mismatch at flat index %d: actual=%d expected=%d\n",
                name, i, actual[i], expected[i]);
            pari_err(e_MISC, "mild certificate: expected matrix mismatch");
        }
    pari_printf("  %s entry-by-entry comparison = OK\n", name);
}

static void
rel_left_multiply3(
    const int left[3][3], int right[REL_ROWS][REL_WORDS],
    int output[REL_ROWS][REL_WORDS])
{
    for (int i = 0; i < REL_ROWS; ++i)
        for (int word = 0; word < REL_WORDS; ++word)
        {
            int value = 0;
            for (int j = 0; j < REL_ROWS; ++j)
                value += left[i][j] * right[j][word];
            output[i][word] = rel_mod(value);
        }
}

static void
rel_substituted_nested_bracket(
    const int degree_one[3][3], int x, int y, int z,
    int output[REL_WORDS])
{
    memset(output, 0, sizeof(int) * REL_WORDS);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
            {
                int coefficient =
                    degree_one[i][x] * degree_one[j][y] * degree_one[k][z]
                    - degree_one[i][y] * degree_one[j][x] * degree_one[k][z]
                    - degree_one[i][z] * degree_one[j][x] * degree_one[k][y]
                    + degree_one[i][z] * degree_one[j][y] * degree_one[k][x];
                output[rel_encode_word(i, j, k)] = rel_mod(coefficient);
            }
}

static void
rel_verify_phi_through_lie(
    const int degree_one[3][3],
    int expected[REL_ROWS][REL_WORDS])
{
    const int triples[REL_LIE_DIM][3] = {
        {0,1,0}, {0,1,1}, {0,1,2}, {0,2,0},
        {0,2,1}, {0,2,2}, {1,2,1}, {1,2,2}
    };
    const int Q[REL_ROWS][REL_LIE_DIM] = {
        {0,0,3,0,4,1,1,4},
        {2,1,4,2,4,2,3,0},
        {4,2,3,0,2,0,1,1}
    };
    int transformed_lie[REL_LIE_DIM][REL_WORDS];
    for (int q = 0; q < REL_LIE_DIM; ++q)
        rel_substituted_nested_bracket(
            degree_one, triples[q][0], triples[q][1], triples[q][2],
            transformed_lie[q]);
    int reconstructed[REL_ROWS][REL_WORDS] = {{0}};
    for (int row = 0; row < REL_ROWS; ++row)
        for (int q = 0; q < REL_LIE_DIM; ++q)
            for (int word = 0; word < REL_WORDS; ++word)
                reconstructed[row][word] =
                    rel_mod(
                        reconstructed[row][word]
                        + Q[row][q] * transformed_lie[q][word]);
    rel_compare_matrix(
        &reconstructed[0][0], &expected[0][0],
        REL_ROWS, REL_WORDS, "independent Lie T_phi");
}

static void
rel_print_polynomial(const char *name, const int row[REL_WORDS])
{
    pari_printf("  %s = ", name);
    int first = 1;
    for (int word = 0; word < REL_WORDS; ++word)
        if (row[word])
        {
            char text[4];
            rel_word_string(word, text);
            pari_printf(
                "%s%s%s", first ? "" : " + ",
                row[word] == 1 ? "" : pari_sprintf("%d ", row[word]),
                text);
            first = 0;
        }
    if (first) pari_printf("0");
    pari_printf("\n");
}

static int
rel_largest_word(const int row[REL_WORDS])
{
    for (int word = REL_WORDS - 1; word >= 0; --word)
        if (row[word]) return word;
    return -1;
}

static void
rel_direct_certificate_path(
    int T[REL_ROWS][REL_WORDS], const int degree_one[3][3],
    const int U[3][3])
{
    int direct_phi[REL_ROWS][REL_WORDS] = {{0}};
    for (int row = 0; row < REL_ROWS; ++row)
        for (int source = 0; source < REL_WORDS; ++source)
        {
            int source_word[3];
            rel_decode_word(source, source_word);
            for (int x = 0; x < 3; ++x)
                for (int y = 0; y < 3; ++y)
                    for (int z = 0; z < 3; ++z)
                        direct_phi[row][rel_encode_word(x,y,z)] =
                            rel_mod(
                                direct_phi[row][rel_encode_word(x,y,z)]
                                + T[row][source]
                                    * degree_one[x][source_word[0]]
                                    * degree_one[y][source_word[1]]
                                    * degree_one[z][source_word[2]]);
        }
    int direct_result[REL_ROWS][REL_WORDS];
    rel_left_multiply3(U, direct_phi, direct_result);
    int leaders[3];
    for (int row = 0; row < 3; ++row)
        leaders[row] = rel_largest_word(direct_result[row]);
    if (leaders[0] != rel_encode_word(2,2,0)
        || leaders[1] != rel_encode_word(2,1,1)
        || leaders[2] != rel_encode_word(2,1,0)
        || !rel_cf_degree3(leaders, 3))
        pari_err(e_MISC, "mild certificate: fresh end-to-end check failed");
    pari_printf(
        "  fresh end-to-end leaders = [cca,cbb,cba], CF = YES\n");
}

static int
rel_contains_forbidden(const int *word, int length)
{
    const int forbidden[3][3] = {
        {2,2,0}, {2,1,1}, {2,1,0}
    };
    for (int start = 0; start + 3 <= length; ++start)
        for (int q = 0; q < 3; ++q)
            if (word[start] == forbidden[q][0]
                && word[start + 1] == forbidden[q][1]
                && word[start + 2] == forbidden[q][2])
                return 1;
    return 0;
}

static void
rel_hilbert_sanity(void)
{
    long expected[13] = {0}, normal[13] = {0};
    expected[0] = 1;
    for (int degree = 1; degree <= 12; ++degree)
        expected[degree] =
            3 * expected[degree - 1]
            - (degree >= 3 ? 3 * expected[degree - 3] : 0);
    normal[0] = 1;
    for (int degree = 1; degree <= 12; ++degree)
    {
        long count = 1;
        for (int i = 0; i < degree; ++i) count *= 3;
        for (long code = 0; code < count; ++code)
        {
            long current = code;
            int word[12];
            for (int i = degree - 1; i >= 0; --i)
            {
                word[i] = current % 3;
                current /= 3;
            }
            if (!rel_contains_forbidden(word, degree)) ++normal[degree];
        }
        if (normal[degree] != expected[degree])
            pari_err(e_MISC, "mild certificate: Hilbert sanity check failed");
    }
    pari_printf("  Hilbert coefficients degree 0..12 = [");
    for (int i = 0; i <= 12; ++i)
        pari_printf("%s%ld", i ? "," : "", normal[i]);
    pari_printf("] = coefficients of 1/(1-3t+3t^3)\n");
}

void
my_run_mild_certificate_fixture(GEN fixture, GEN p)
{
    if (!equaliu(p, REL_P) || typ(fixture) != t_MAT
        || glength(fixture) != REL_WORDS
        || lg(gel(fixture, 1)) != REL_ROWS + 1)
        pari_err(e_MISC, "mild certificate: fixture must be 3 x 27 over F_5");
    rel_current_p = REL_P;
    int T[REL_ROWS][REL_WORDS];
    for (int column = 0; column < REL_WORDS; ++column)
        for (int row = 0; row < REL_ROWS; ++row)
            T[row][column] =
                (int)umodiu(gcoeff(fixture, row + 1, column + 1), REL_P);
    if (rel_matrix_rank(REL_ROWS, REL_WORDS, T) != 3)
        pari_err(e_MISC, "mild certificate: fixture rank changed");

    /*
     * Column j is the coefficient vector of phi(generator j) in (a,b,c).
     * These entries are constructed directly from:
     * phi(a)=c, phi(b)=b+c, phi(c)=a+3b.
     */
    const int M_phi[3][3] = {
        {0,0,1},
        {0,1,3},
        {1,1,0}
    };
    const int expected_M_inverse[3][3] = {
        {3,4,1},
        {2,1,0},
        {1,0,0}
    };
    int M_inverse[3][3];
    rel_inverse_square(3, &M_phi[0][0], &M_inverse[0][0]);
    pari_printf("MASSEY_CERTIFICATE_TEST\n");
    pari_printf(
        "  column-image convention M_phi = [0,0,1;0,1,3;1,1,0]\n");
    pari_printf("  det(M_phi) = %d\n", rel_det3(M_phi));
    if (rel_det3(M_phi) != 4)
        pari_err(e_MISC, "mild certificate: det(M_phi) is not 4");
    rel_compare_matrix(
        &M_inverse[0][0], &expected_M_inverse[0][0],
        3, 3, "M_phi inverse");
    rel_verify_inverse(
        3, &M_phi[0][0], &M_inverse[0][0], "M_phi");

    int cubic[REL_WORDS][REL_WORDS];
    int cubic_inverse[REL_WORDS][REL_WORDS];
    rel_build_cubic_substitution(M_phi, cubic);
    rel_build_cubic_substitution(
        (const int (*)[3])M_inverse, cubic_inverse);
    int computed_cubic_inverse[REL_WORDS][REL_WORDS];
    rel_inverse_square(
        REL_WORDS, &cubic[0][0], &computed_cubic_inverse[0][0]);
    if (rel_rank27(cubic) != 27)
        pari_err(e_MISC, "mild certificate: cubic substitution rank is not 27");
    rel_compare_matrix(
        &computed_cubic_inverse[0][0], &cubic_inverse[0][0],
        REL_WORDS, REL_WORDS, "degree-3 inverse operator");
    rel_verify_inverse(
        REL_WORDS, &cubic[0][0], &cubic_inverse[0][0],
        "degree-3 phi");
    pari_printf("  rank(degree-3 phi) = 27\n");

    int T_phi[REL_ROWS][REL_WORDS];
    rel_apply_cubic_substitution(T, cubic, T_phi);
    const int expected_T_phi[REL_ROWS][REL_WORDS] = {
        {0,4,0,2,2,2,0,2,0,4,1,1,2,0,4,2,2,0,0,1,0,2,4,0,0,0,0},
        {0,0,2,0,2,2,1,2,1,0,1,1,2,0,2,2,1,0,2,1,3,2,2,0,1,0,0},
        {0,1,1,3,1,0,3,3,2,1,3,2,1,0,4,3,2,0,1,2,1,0,4,0,2,0,0}
    };
    rel_compare_matrix(
        &T_phi[0][0], &expected_T_phi[0][0],
        REL_ROWS, REL_WORDS, "directly computed T_phi");
    rel_print_int_matrix(
        "  T_phi = ", REL_ROWS, REL_WORDS,
        &T_phi[0][0], REL_WORDS);
    if (rel_matrix_rank(REL_ROWS, REL_WORDS, T_phi) != 3)
        pari_err(e_MISC, "mild certificate: rank(T_phi) is not 3");
    int recovered_T[REL_ROWS][REL_WORDS];
    rel_apply_cubic_substitution(T_phi, cubic_inverse, recovered_T);
    rel_compare_matrix(
        &recovered_T[0][0], &T[0][0],
        REL_ROWS, REL_WORDS, "phi inverse recovery of T");
    rel_verify_phi_through_lie(M_phi, T_phi);

    const int U[3][3] = {
        {2,3,4},
        {4,1,2},
        {0,3,1}
    };
    int U_inverse[3][3];
    rel_inverse_square(3, &U[0][0], &U_inverse[0][0]);
    pari_printf("  U = [2,3,4;4,1,2;0,3,1]\n");
    pari_printf("  det(U) = %d\n", rel_det3(U));
    if (rel_det3(U) != 1)
        pari_err(e_MISC, "mild certificate: det(U) is not 1");
    rel_verify_inverse(3, &U[0][0], &U_inverse[0][0], "U");

    int R_cert[REL_ROWS][REL_WORDS];
    rel_left_multiply3(U, T_phi, R_cert);
    const int expected_R[REL_ROWS][REL_WORDS] = {
        {0,2,0,1,4,0,0,2,1,2,2,3,4,0,0,2,0,0,0,3,3,0,0,0,1,0,0},
        {0,3,4,4,2,0,2,1,0,3,1,4,2,0,1,1,3,0,4,4,0,0,1,0,0,0,0},
        {0,1,2,3,2,1,1,4,0,1,1,0,2,0,0,4,0,0,2,0,0,1,0,0,0,0,0}
    };
    rel_compare_matrix(
        &R_cert[0][0], &expected_R[0][0],
        REL_ROWS, REL_WORDS, "computed R_cert");
    rel_print_int_matrix(
        "  R_cert = ", REL_ROWS, REL_WORDS,
        &R_cert[0][0], REL_WORDS);
    int recovered_phi[REL_ROWS][REL_WORDS];
    rel_left_multiply3(
        (const int (*)[3])U_inverse, R_cert, recovered_phi);
    rel_compare_matrix(
        &recovered_phi[0][0], &T_phi[0][0],
        REL_ROWS, REL_WORDS, "U inverse recovery of T_phi");

    rel_print_polynomial("rho1", R_cert[0]);
    rel_print_polynomial("rho2", R_cert[1]);
    rel_print_polynomial("rho3", R_cert[2]);

    const int expected_leaders[3] = {
        24, 22, 21
    };
    int leaders[3];
    pari_printf("  leading means lexicographically largest for a<b<c\n");
    for (int row = 0; row < 3; ++row)
    {
        leaders[row] = rel_largest_word(R_cert[row]);
        char leader[4];
        rel_word_string(leaders[row], leader);
        pari_printf(
            "  LT(rho%d) = %s coefficient=%d; larger coefficients:",
            row + 1, leader, R_cert[row][leaders[row]]);
        for (int word = leaders[row] + 1; word < REL_WORDS; ++word)
        {
            char text[4];
            rel_word_string(word, text);
            pari_printf(" %s=%d", text, R_cert[row][word]);
        }
        pari_printf("\n");
        if (leaders[row] != expected_leaders[row]
            || R_cert[row][leaders[row]] != 1)
            pari_err(e_MISC, "mild certificate: leading word mismatch");
    }
    if (!rel_cf_degree3(leaders, 3))
        pari_err(e_MISC, "mild certificate: leaders are not CF");
    pari_printf("  overlap certificate (prefix of first vs suffix of second):\n");
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            int left[3], right[3];
            char wi[4], wj[4];
            rel_decode_word(leaders[i], left);
            rel_decode_word(leaders[j], right);
            rel_word_string(leaders[i], wi);
            rel_word_string(leaders[j], wj);
            for (int length = 1; length <= 2; ++length)
            {
                int equal = 1;
                for (int q = 0; q < length; ++q)
                    if (left[q] != right[3 - length + q])
                        equal = 0;
                pari_printf(
                    "    %s prefix(%d) vs %s suffix(%d): %s\n",
                    wi, length, wj, length,
                    equal ? "EQUAL (FAIL)" : "different");
                if (equal)
                    pari_err(e_MISC, "mild certificate: overlap found");
            }
        }
    pari_printf(
        "  factor condition: equal length and pairwise distinct = OK\n");
    pari_printf("  combinatorially free = YES\n");

    rel_direct_certificate_path(T, M_phi, U);
    rel_hilbert_sanity();

    pari_printf(
        "  Anick criterion: homogeneous cubic relations with "
        "combinatorially free leading monomials are strongly free = APPLIES\n");
    pari_printf(
        "  U invariance: U is invertible, so both bases generate the same "
        "homogeneous ideal with the same three minimal degree-3 generators\n");
    pari_printf(
        "  phi invariance: phi is a graded algebra automorphism; it maps "
        "the generated ideals and quotients isomorphically, preserving "
        "Hilbert series and strong freeness\n");
    pari_printf(
        "MASSEY_CERTIFICATE_TEST completed without number-field arithmetic\n");
    pari_printf(
        "VERIFIED CERTIFICATE: leaders [cca,cbb,cba]; "
        "R3 strongly free; tower group mild; cd(G)=2\n");
}

static int
rel_vector_weight(const int vector[3])
{
    return (vector[0] != 0) + (vector[1] != 0) + (vector[2] != 0);
}

static int
rel_vector_compare(const void *left, const void *right)
{
    const int *x = left, *y = right;
    int wx = rel_vector_weight(x), wy = rel_vector_weight(y);
    if (wx != wy) return wx - wy;
    for (int i = 0; i < 3; ++i)
        if (x[i] != y[i]) return x[i] - y[i];
    return 0;
}

int
my_anick_words_combinatorially_free(GEN words)
{
    if ((typ(words) != t_VEC && typ(words) != t_COL)
        || glength(words) != REL_ROWS)
        pari_err_TYPE("Anick leading words", words);
    int encoded[REL_ROWS];
    for (int i = 0; i < REL_ROWS; ++i)
    {
        if (typ(gel(words, i + 1)) != t_INT)
            pari_err_TYPE("Anick encoded word", gel(words, i + 1));
        encoded[i] = itos(gel(words, i + 1));
        if (encoded[i] < 0 || encoded[i] >= REL_WORDS)
            pari_err_DOMAIN(
                "Anick encoded word", "word", "outside", stoi(REL_WORDS),
                gel(words, i + 1));
    }
    return rel_cf_degree3(encoded, REL_ROWS);
}

GEN
my_find_strongly_free_witness(
    GEN fixture, GEN p, long candidate_limit)
{
    pari_sp av = avma;
    if (typ(p) != t_INT || !uisprime(itou(p)) || !mpodd(p))
        pari_err_TYPE("strong-freeness search [odd prime]", p);
    long prime = itos(p);
    if (prime > 31)
        pari_err(e_MISC, "strong-freeness search supports primes at most 31");
    if (typ(fixture) != t_MAT || glength(fixture) != REL_WORDS
        || lg(gel(fixture, 1)) != REL_ROWS + 1)
        pari_err(e_MISC, "strong-freeness search requires a 3 x 27 matrix");
    if (candidate_limit == 0)
        pari_err(e_MISC, "strong-freeness candidate limit must be positive");

    int previous_prime = rel_current_p;
    rel_current_p = (int)prime;
    int T[REL_ROWS][REL_WORDS];
    for (int column = 0; column < REL_WORDS; ++column)
        for (int row = 0; row < REL_ROWS; ++row)
            T[row][column] =
                (int)umodiu(gcoeff(fixture, row + 1, column + 1), prime);
    if (rel_matrix_rank(REL_ROWS, REL_WORDS, T) != REL_ROWS)
    {
        rel_current_p = previous_prime;
        return gerepilecopy(av, cgetg(1, t_VEC));
    }

    long vector_count = prime * prime * prime - 1;
    int (*vectors)[3] =
        malloc((size_t)vector_count * sizeof(*vectors));
    if (!vectors) pari_err(e_MEM);
    long at = 0;
    for (int x = 0; x < prime; ++x)
        for (int y = 0; y < prime; ++y)
            for (int z = 0; z < prime; ++z)
                if (x || y || z)
                {
                    vectors[at][0] = x;
                    vectors[at][1] = y;
                    vectors[at][2] = z;
                    ++at;
                }
    qsort(vectors, (size_t)vector_count, sizeof(*vectors), rel_vector_compare);

    long effective_limit = candidate_limit;
    if (candidate_limit < 0)
    {
        long q = prime;
        effective_limit =
            (q * q * q - 1)
            * (q * q * q - q)
            * (q * q * q - q * q);
    }
    long candidates = 0;
    int weights[3] = {1, 1, 1};
    for (long a = 0; a < vector_count; ++a)
        for (long b = 0; b < vector_count; ++b)
            for (long c = 0; c < vector_count; ++c)
            {
                int M[3][3] = {
                    {vectors[a][0], vectors[b][0], vectors[c][0]},
                    {vectors[a][1], vectors[b][1], vectors[c][1]},
                    {vectors[a][2], vectors[b][2], vectors[c][2]}
                };
                if (!rel_det3(M)) continue;
                if (++candidates > effective_limit)
                {
                    free(vectors);
                    rel_current_p = previous_prime;
                    return gerepilecopy(av, cgetg(1, t_VEC));
                }

                int cubic[REL_WORDS][REL_WORDS];
                int transformed[REL_ROWS][REL_WORDS];
                rel_build_cubic_substitution(M, cubic);
                rel_apply_cubic_substitution(T, cubic, transformed);
                for (int permutation = 0; permutation < 6; ++permutation)
                {
                    int order[REL_WORDS], reduced[REL_ROWS][REL_WORDS];
                    int operations[REL_ROWS][REL_ROWS], pivots[REL_ROWS];
                    rel_make_order(
                        rel_permutations[permutation], weights, order);
                    rel_ordered_rref(
                        transformed, order, reduced, operations, pivots);
                    if (!rel_cf_degree3(pivots, REL_ROWS)) continue;

                    GEN M_gen = zeromatcopy(3, 3);
                    GEN U_gen = zeromatcopy(3, 3);
                    GEN R_gen = zeromatcopy(REL_ROWS, REL_WORDS);
                    GEN leaders = cgetg(REL_ROWS + 1, t_VEC);
                    GEN variable_order = cgetg(4, t_VEC);
                    for (int i = 0; i < 3; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            gcoeff(M_gen, i + 1, j + 1) = stoi(M[i][j]);
                            gcoeff(U_gen, i + 1, j + 1) =
                                stoi(operations[i][j]);
                        }
                        char word[4];
                        rel_word_string(pivots[i], word);
                        gel(leaders, i + 1) = strtoGENstr(word);
                        char letter[2] = {
                            rel_letters[rel_permutations[permutation][i]], '\0'
                        };
                        gel(variable_order, i + 1) = strtoGENstr(letter);
                    }
                    for (int row = 0; row < REL_ROWS; ++row)
                        for (int column = 0; column < REL_WORDS; ++column)
                            gcoeff(R_gen, row + 1, column + 1) =
                                stoi(reduced[row][column]);
                    free(vectors);
                    rel_current_p = previous_prime;
                    GEN witness = cgetg(7, t_VEC);
                    gel(witness, 1) = M_gen;
                    gel(witness, 2) = U_gen;
                    gel(witness, 3) = R_gen;
                    gel(witness, 4) = leaders;
                    gel(witness, 5) = variable_order;
                    gel(witness, 6) = gen_1;
                    return gerepilecopy(av, witness);
                }
            }

    free(vectors);
    rel_current_p = previous_prime;
    return gerepilecopy(av, cgetg(1, t_VEC));
}
