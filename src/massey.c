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


#include <pthread.h>
#include <pari/pari.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// Debug level
#define MY_DEBUGLEVEL 0

// Debug printing function
#define DEBUG_PRINT(level, ...) \
    do { if (MY_DEBUGLEVEL >= (level)) pari_printf(__VA_ARGS__); } while (0)

#include "../headers/misc_functions.h"
#include "../headers/artin_symbol.h"
#include "../headers/ext_and_aut.h"
#include "../headers/find_cup_matrix.h"
#include "../headers/secondary_norm.h"
#include "../headers/massey_tensor.h"
#include "../headers/relation_degree3.h"

// Function prototype for parallel computation
GEN compute_my_relations(long i, GEN args);

static GEN
my_evaluate_character_on_operator(GEN p, GEN character, GEN D)
{
    pari_sp av = avma;
    long inputs = lg(D) - 1;
    GEN result = zerocol(inputs);
    for (long j = 1; j <= inputs; ++j)
        gel(result, j) =
            FpV_dotproduct(character, gel(D, j), p);
    return gerepilecopy(av, result);
}

static void
my_validate_p_rank_helper(void)
{
    GEN p5 = stoi(5);
    GEN p3 = stoi(3);
    GEN cases[] = {
        mkvec2(stoi(10), stoi(5)),
        mkvec2(stoi(130), stoi(10)),
        mkvec3(stoi(25), stoi(10), stoi(7)),
        mkvec3(stoi(125), stoi(25), stoi(5)),
        mkvec3(stoi(2), stoi(3), stoi(7)),
        mkvec3(stoi(9), stoi(3), stoi(2))
    };
    GEN primes[] = { p5, p5, p5, p5, p5, p3 };
    long expected[] = { 2, 2, 2, 3, 0, 2 };

    for (long i = 0; i < 6; ++i)
    {
        long actual = my_p_rank_from_cyc(cases[i], primes[i]);
        if (actual != expected[i])
            pari_err(e_MISC, "synthetic p-rank regression failed");
        pari_printf(
            "  p-rank(%Ps, p=%Ps) = %ld, OK\n",
            cases[i], primes[i], actual);
    }
}

static double
my_elapsed_seconds(struct timespec started, struct timespec finished)
{
    return (double)(finished.tv_sec - started.tv_sec)
        + (double)(finished.tv_nsec - started.tv_nsec) / 1e9;
}

static GEN
my_audit_character_kernel(GEN K, GEN p, GEN t)
{
    GEN cyc = bnf_get_cyc(K);
    GEN chi = cgetg(lg(cyc), t_VEC);
    long coordinate = 1;
    for (long i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
            gel(chi, i) =
                modii(
                    mulii(
                        diviiexact(gel(cyc, i), p),
                        gel(t, coordinate++)),
                    gel(cyc, i));
        else
            gel(chi, i) = gen_0;
    return charker(cyc, chi);
}

static void
my_run_arithmetic_audit(
    GEN K, GEN p, GEN discriminant, GEN Ja_vect, GEN D_prime_vect)
{
    const char *sample_names[] = { "a", "b", "c", "a+b", "a+c", "b+c" };
    const char *coordinate_names[] = { "a", "b", "c" };
    GEN characters[] = {
        mkcol3(gen_1, gen_0, gen_0),
        mkcol3(gen_0, gen_1, gen_0),
        mkcol3(gen_0, gen_0, gen_1),
        mkcol3(gen_1, gen_1, gen_0),
        mkcol3(gen_1, gen_0, gen_1),
        mkcol3(gen_0, gen_1, gen_1)
    };

    if (!equalii(discriminant, stoi(-11203620)) || !equaliu(p, 5))
        pari_err(e_MISC, "MASSEY_ARITHMETIC_AUDIT requires the fixed p=5 field");
    if (!gequal(bnf_get_cyc(K), mkvec3(stoi(10), stoi(10), stoi(10))))
        pari_err(e_MISC, "audit class-group invariants are not [10,10,10]");
    if (glength(Ja_vect) != 3)
        pari_err(e_MISC, "audit Ja_vect does not have three entries");

    pari_printf("\nMASSEY_ARITHMETIC_AUDIT base and Ja certificate\n");
    pari_printf("  p=5; disc(K)=%Ps; Cl(K)=%Ps\n", discriminant, bnf_get_cyc(K));
    GEN torsion_coordinates = cgetg(4, t_MAT);
    for (long j = 1; j <= 3; ++j)
    {
        GEN a_prime = gmael(Ja_vect, j, 1);
        GEN J = gmael(Ja_vect, j, 2);
        GEN relation =
            idealmul(
                K, idealhnf0(K, a_prime, NULL),
                idealpow(K, J, p));
        if (!gequal(
                idealhnf0(K, relation, NULL),
                idealhnf0(K, gen_1, NULL)))
            pari_err(e_MISC, "audit div(a')+5J is not zero");
        GEN class_exp = bnfisprincipal0(K, J, 0);
        GEN killed = bnfisprincipal0(K, idealpow(K, J, p), 0);
        if (!ZV_equal0(killed))
            pari_err(e_MISC, "audit Ja class is not killed by 5");
        GEN coordinates = cgetg(4, t_COL);
        for (long i = 1; i <= 3; ++i)
        {
            GEN residue = modii(gel(class_exp, i), stoi(10));
            if (signe(modii(residue, stoi(2))))
                pari_err(e_MISC, "audit Ja class is not in Cl(K)[5]");
            gel(coordinates, i) = modii(diviiexact(residue, stoi(2)), p);
        }
        gel(torsion_coordinates, j) = coordinates;
        pari_printf(
            "  e_%ld=(a'=%Ps, J=%Ps): div(a')+5J=0 PASS; "
            "[J]=%Ps; Cl(K)[5] coordinates=%Ps; killed-by-5 PASS\n",
            j, a_prime, J, gtovec(class_exp), gtovec(coordinates));
    }
    if (FpM_rank(torsion_coordinates, p) != 3)
        pari_err(e_MISC, "audit Ja classes are not F_5-independent");
    pari_printf(
        "  Cl(K)[5] coordinate matrix=%Ps; rank=3; Ja_vect basis PASS\n",
        torsion_coordinates);

    GEN family = my_secondary_norm_basis_family(K, p, Ja_vect, D_prime_vect);
    GEN samples[6];
    for (long i = 0; i < 3; ++i) samples[i] = gmael(family, 2, i + 1);
    for (long i = 0; i < 3; ++i) samples[i + 3] = gmael(family, 3, i + 1);

    GEN expected[] = {
        gp_read_str("[0,0,0;0,3,1;0,3,0]"),
        gp_read_str("[0,1,2;0,0,0;4,2,4]"),
        gp_read_str("[1,2,0;4,0,1;0,0,0]"),
        gp_read_str("[0,3,1;0,2,4;0,1,2]"),
        gp_read_str("[1,4,0;1,4,4;4,1,0]"),
        gp_read_str("[3,1,2;0,3,2;0,2,3]")
    };
    pari_printf("\nMASSEY_ARITHMETIC_AUDIT final matrix comparison\n");
    for (long q = 0; q < 6; ++q)
    {
        if (!gequal(samples[q], expected[q]))
            pari_err(e_MISC, "audited matrix differs from expected fixture");
        pari_printf("  D_(%s)=%Ps EXACT MATCH\n", sample_names[q], samples[q]);
    }

    GEN T = my_triple_massey_word_matrix(p, family);
    pari_printf("\nMASSEY_ARITHMETIC_AUDIT 54 final Q cross-checks\n");
    for (long q = 0; q < 6; ++q)
        for (long y = 0; y < 3; ++y)
        {
            GEN direct =
                my_evaluate_character_on_operator(p, characters[y], samples[q]);
            GEN tensor =
                my_triple_massey_contract(
                    p, T, 3, characters[q], characters[q], characters[y]);
            if (!gequal(direct, tensor))
                pari_err(e_MISC, "audit final tensor cross-check failed");
            for (long j = 1; j <= 3; ++j)
                pari_printf(
                    "  x=%s y=%s e_%ld: direct=%Ps tensor=%Ps PASS\n",
                    sample_names[q], coordinate_names[y], j,
                    gel(direct, j), gel(tensor, j));
        }

    pari_printf("\nMASSEY_ARITHMETIC_AUDIT full scalar arithmetic\n");
    for (long i = 0; i < 3; ++i)
    {
        GEN twice = zerocol(3);
        gel(twice, i + 1) = stoi(2);
        GEN H1 = my_audit_character_kernel(K, p, characters[i]);
        GEN H2 = my_audit_character_kernel(K, p, twice);
        if (!gequal(H1, H2))
            pari_err(e_MISC, "audit scalar multiple changed class-field kernel");
        GEN D_twice =
            my_secondary_norm_operator(K, p, twice, Ja_vect, D_prime_vect);
        GEN four_D = FpM_Fp_mul(samples[i], stoi(4), p);
        if (!gequal(D_twice, four_D))
            pari_err(e_MISC, "audit full arithmetic quadratic scaling failed");
        pari_printf(
            "  H_(%s)=H_(2%s)=%Ps; D_(2%s)=%Ps=4D_(%s) PASS\n",
            coordinate_names[i], coordinate_names[i], H1,
            coordinate_names[i], D_twice, coordinate_names[i]);
    }
    pari_printf("\nARITHMETIC CERTIFICATE VERIFIED\n");
}

static long
my_rank3_word_index(long i, long j, long k)
{
    return ((i - 1) * 3 + (j - 1)) * 3 + k;
}

static void
my_rank3_validate_tensor(GEN p, GEN T)
{
    for (long i = 1; i <= 3; ++i)
        for (long j = 1; j <= 3; ++j)
            for (long k = 1; k <= 3; ++k)
            {
                GEN Mijk = gel(T, my_rank3_word_index(i, j, k));
                GEN Mkji = gel(T, my_rank3_word_index(k, j, i));
                if (!gequal(Mijk, Mkji))
                {
                    pari_printf(
                        "rank-3 outer symmetry failure: "
                        "(%ld,%ld,%ld)=%Ps, (%ld,%ld,%ld)=%Ps\n",
                        i, j, k, gtovec(Mijk),
                        k, j, i, gtovec(Mkji));
                    pari_err(e_MISC, "rank-3 outer-symmetry regression failed");
                }

                GEN Mjki = gel(T, my_rank3_word_index(j, k, i));
                GEN Mkij = gel(T, my_rank3_word_index(k, i, j));
                GEN cyclic = FpV_add(FpV_add(Mijk, Mjki, p), Mkij, p);
                if (!gequal0(cyclic))
                {
                    pari_printf(
                        "rank-3 cyclic failure at (%ld,%ld,%ld): "
                        "%Ps + %Ps + %Ps = %Ps\n",
                        i, j, k, gtovec(Mijk), gtovec(Mjki),
                        gtovec(Mkij), gtovec(cyclic));
                    pari_err(e_MISC, "rank-3 cyclic-shuffle regression failed");
                }

                if (i == j && j == k && !gequal0(Mijk))
                {
                    pari_printf(
                        "rank-3 diagonal failure at (%ld,%ld,%ld): %Ps\n",
                        i, j, k, gtovec(Mijk));
                    pari_err(e_MISC, "rank-3 diagonal regression failed");
                }
            }
    my_validate_triple_massey_identities(p, T, 3);
}

static void
my_run_rank3_test(
    GEN K, GEN p, GEN discriminant, GEN Kcyc,
    GEN Ja_vect, GEN D_prime_vect,
    double bnf_time, double ja_time,
    struct timespec total_started)
{
    const char *names[] = { "a", "b", "c" };
    const char *sample_names[] = { "D_a", "D_b", "D_c", "D_ab", "D_ac", "D_bc" };
    GEN characters[] = {
        mkcol3(gen_1, gen_0, gen_0),
        mkcol3(gen_0, gen_1, gen_0),
        mkcol3(gen_0, gen_0, gen_1),
        mkcol3(gen_1, gen_1, gen_0),
        mkcol3(gen_1, gen_0, gen_1),
        mkcol3(gen_0, gen_1, gen_1)
    };
    long p_rank = my_p_class_rank(K, p);

    pari_printf("\nMASSEY_RANK3_TEST initialization\n");
    pari_printf("  discriminant = %Ps\n", discriminant);
    pari_printf("  PARI class-group invariants = %Ps\n", Kcyc);
    pari_printf("  p-rank = %ld\n", p_rank);
    pari_printf("  Ja_vect length = %ld\n", glength(Ja_vect));
    pari_printf("  D_prime_vect length = %ld\n", glength(D_prime_vect));
    pari_printf("  base BNF time = %.3f s\n", bnf_time);
    pari_printf("  Ja_vect preparation time = %.3f s\n", ja_time);
    pari_printf(
        "  basis convention: a=[1,0,0], b=[0,1,0], c=[0,0,1]\n");
    pari_printf(
        "  sample order: [1,0,0], [0,1,0], [0,0,1], "
        "[1,1,0], [1,0,1], [0,1,1]\n");
    fflush(stdout);

    if (!equalii(discriminant, stoi(-11203620)))
        pari_err(e_MISC, "MASSEY_RANK3_TEST requires discriminant -11203620");
    if (!equaliu(p, 5))
        pari_err(e_MISC, "MASSEY_RANK3_TEST requires p=5");
    if (p_rank != 3)
        pari_err(e_MISC, "rank-3 field does not have 5-class rank 3");
    if (glength(Ja_vect) != 3)
        pari_err(e_MISC, "rank-3 Ja_vect does not have length 3");

    struct timespec algebra_started, algebra_finished;
    GEN family =
        my_secondary_norm_basis_family(K, p, Ja_vect, D_prime_vect);
    if (itos(gel(family, 1)) != 3)
        pari_err(e_MISC, "rank-3 basis family has the wrong dimension");
    clock_gettime(CLOCK_MONOTONIC, &algebra_started);

    GEN samples[6];
    for (long i = 0; i < 3; ++i)
        samples[i] = gel(gel(family, 2), i + 1);
    for (long i = 0; i < 3; ++i)
        samples[i + 3] = gel(gel(family, 3), i + 1);

    pari_printf("\nMASSEY_RANK3_TEST arithmetic family\n");
    for (long q = 0; q < 6; ++q)
    {
        GEN Dq = samples[q];
        if (typ(Dq) != t_MAT || glength(Dq) != 3
            || lg(gel(Dq, 1)) != 4)
            pari_err(e_MISC, "rank-3 secondary-norm matrix is not 3 x 3");
        for (long col = 1; col <= 3; ++col)
            for (long row = 1; row <= 3; ++row)
            {
                GEN entry = gcoeff(Dq, row, col);
                if (typ(entry) != t_INT || signe(entry) < 0
                    || cmpii(entry, p) >= 0)
                    pari_err(e_MISC,
                             "rank-3 secondary-norm entry is not in F_5");
            }
        GEN reconstructed =
            my_reconstruct_secondary_norm(p, family, characters[q]);
        if (!gequal(reconstructed, Dq))
            pari_err(e_MISC, "rank-3 sampled-character reconstruction failed");
        pari_printf(
            "  %s = D_%Ps = %Ps\n",
            sample_names[q], gtovec(characters[q]), Dq);
    }

    GEN zero = zerocol(3);
    if (!gequal0(my_reconstruct_secondary_norm(p, family, zero)))
        pari_err(e_MISC, "rank-3 zero-character reconstruction failed");
    for (long i = 0; i < 3; ++i)
    {
        GEN twice = zerocol(3);
        gel(twice, i + 1) = stoi(2);
        GEN reconstructed =
            my_reconstruct_secondary_norm(p, family, twice);
        GEN scaled = FpM_Fp_mul(samples[i], stoi(4), p);
        if (!gequal(reconstructed, scaled))
            pari_err(e_MISC, "rank-3 quadratic scaling failed");
    }

    pari_printf("\nMASSEY_RANK3_TEST basis polarizations\n");
    for (long i = 1; i <= 3; ++i)
        for (long k = i; k <= 3; ++k)
        {
            GEN delta =
                my_secondary_norm_delta_basis(p, family, i, k);
            GEN reverse =
                my_secondary_norm_delta_basis(p, family, k, i);
            if (!gequal(delta, reverse))
                pari_err(e_MISC, "rank-3 DeltaD symmetry failed");
            pari_printf(
                "  DeltaD(%s,%s) = %Ps\n",
                names[i - 1], names[k - 1], delta);
        }

    GEN T = my_triple_massey_word_matrix(p, family);
    if (typ(T) != t_MAT || glength(T) != 27 || lg(gel(T, 1)) != 4)
        pari_err(e_MISC, "rank-3 tensor matrix is not 3 x 27");
    my_rank3_validate_tensor(p, T);
    long tensor_rank = FpM_rank(T, p);

    pari_printf("\nMASSEY_RANK3_TEST complete word matrix\n");
    pari_printf("  T = %Ps\n", T);
    pari_printf("  labeled columns:\n");
    for (long i = 1; i <= 3; ++i)
        for (long j = 1; j <= 3; ++j)
            for (long k = 1; k <= 3; ++k)
                pari_printf(
                    "    %2ld -> %s%s%s\n",
                    my_rank3_word_index(i, j, k),
                    names[i - 1], names[j - 1], names[k - 1]);

    pari_printf("  nonzero word evaluations:\n");
    for (long i = 1; i <= 3; ++i)
        for (long j = 1; j <= 3; ++j)
            for (long k = 1; k <= 3; ++k)
            {
                GEN value = gel(T, my_rank3_word_index(i, j, k));
                if (!gequal0(value))
                    pari_printf(
                        "    M(%s,%s,%s) = %Ps\n",
                        names[i - 1], names[j - 1], names[k - 1],
                        gtovec(value));
            }

    pari_printf("\nMASSEY_RANK3_TEST Q checks\n");
    for (long x = 0; x < 3; ++x)
        for (long y = 0; y < 3; ++y)
        {
            GEN direct =
                my_evaluate_character_on_operator(
                    p, characters[y], samples[x]);
            GEN tensor =
                my_triple_massey_contract(
                    p, T, 3, characters[x], characters[x], characters[y]);
            if (!gequal(direct, tensor))
                pari_err(e_MISC, "rank-3 Q(x,y) regression failed");
            pari_printf(
                "  Q(%s,%s): direct=%Ps tensor=%Ps, OK\n",
                names[x], names[y], gtovec(direct), gtovec(tensor));
        }

    GEN x = characters[3];
    GEN y = characters[5];
    GEN z = characters[4];
    GEN x_plus_z = FpV_add(x, z, p);
    GEN D_x = my_reconstruct_secondary_norm(p, family, x);
    GEN D_z = my_reconstruct_secondary_norm(p, family, z);
    GEN D_x_plus_z =
        my_reconstruct_secondary_norm(p, family, x_plus_z);
    GEN polarization =
        FpM_sub(FpM_add(D_x, D_z, p), D_x_plus_z, p);
    GEN direct =
        my_evaluate_character_on_operator(p, y, polarization);
    GEN tensor =
        my_triple_massey_contract(p, T, 3, x, y, z);
    if (!gequal(direct, tensor))
        pari_err(e_MISC, "rank-3 non-basis contraction regression failed");

    clock_gettime(CLOCK_MONOTONIC, &algebra_finished);
    struct timespec total_finished;
    clock_gettime(CLOCK_MONOTONIC, &total_finished);
    pari_printf("\nMASSEY_RANK3_TEST summary\n");
    pari_printf("  rank(T) = %ld\n", tensor_rank);
    pari_printf("  row-span dimension = %ld\n", tensor_rank);
    pari_printf(
        "  three rows are %s\n",
        tensor_rank == 3 ? "independent" : "dependent");
    pari_printf("  tensor identities = OK\n");
    pari_printf(
        "  non-basis M(a+b,b+c,a+c): direct=%Ps tensor=%Ps, OK\n",
        gtovec(direct), gtovec(tensor));
    pari_printf(
        "  tensor/algebra time = %.3f s\n",
        my_elapsed_seconds(algebra_started, algebra_finished));
    pari_printf(
        "  total rank-3 test time = %.3f s\n",
        my_elapsed_seconds(total_started, total_finished));
    pari_printf("MASSEY_RANK3_TEST completed successfully\n");
    fflush(stdout);
}

int
main (int argc, char *argv[])	  
{
    const char *certificate_test = getenv("MASSEY_CERTIFICATE_TEST");
    if (certificate_test && strcmp(certificate_test, "1") == 0)
    {
        pari_init(4000000, 500000);
        GEN fixture =
            gp_read_str(
                "[0,0,0,0,0,3,0,4,1,0,0,3,0,0,4,4,2,4,0,3,3,3,4,2,1,4,0;"
                "0,3,3,4,1,4,4,4,2,3,3,2,1,0,2,4,1,0,3,2,1,4,2,0,2,0,0;"
                "0,1,0,3,2,3,0,2,0,1,1,0,2,0,4,2,2,1,0,0,0,3,4,3,0,1,0]");
        my_run_mild_certificate_fixture(fixture, stoi(5));
        pari_close();
        return 0;
    }

    const char *relation_test = getenv("MASSEY_RELATION_TEST");
    if (relation_test && strcmp(relation_test, "1") == 0)
    {
        pari_init(4000000, 500000);
        GEN fixture =
            gp_read_str(
                "[0,0,0,0,0,3,0,4,1,0,0,3,0,0,4,4,2,4,0,3,3,3,4,2,1,4,0;"
                "0,3,3,4,1,4,4,4,2,3,3,2,1,0,2,4,1,0,3,2,1,4,2,0,2,0,0;"
                "0,1,0,3,2,3,0,2,0,1,1,0,2,0,4,2,2,1,0,0,0,3,4,3,0,1,0]");
        my_run_relation_degree3_fixture(fixture, stoi(5));
        pari_close();
        return 0;
    }

    if (argc != 3) {
        fprintf(stderr,
                "Usage: %s <prime p> <defining polynomial>\n"
                "Example: %s 3 \"x^2+4027\"\n",
                argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    printf(ANSI_COLOR_YELLOW "\n---------------------------------------------------------------------------------------------------------\nStarting program: Massey products\n---------------------------------------------------------------------------------------------------------\n\n" ANSI_COLOR_RESET);

    // Start timer (actual time)
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    
    // Start timer (CPU time)
    clock_t start = clock();

    int min, sec, msec;
    
    //--------------------------------------------------
    // Initialize PARI/GP
    //--------------------------------------------------
    // pari_init(1L<<30,500000);
    entree ep = {"_worker",0,(void*)compute_my_relations,20,"LG",""};
    pari_init_opts(1L<<30,1048576, INIT_JMPm|INIT_SIGm|INIT_DFTm|INIT_noIMTm);
    pari_add_function(&ep); /* add Cworker function to gp */
    pari_mt_init(); /* ... THEN initialize parallelism */
    paristack_setsize(1L<<30, 1L<<33);
    sd_threadsizemax("2147483648", 0);
    setalldebug(0);
    //--------------------------------------------------
    
    int p_int, p_rk, r_rk;
    GEN p, K, f, Kcyc, p_ClFld_pol, J_vect, Ja_vect, D, D_prime_vect;

    // Read the prime number p from arguments
    p = gp_read_str(argv[1]);
    p_int = atoi(argv[1]);
    
    // Read the defining polynomial for K
    f = gp_read_str(argv[2]);
    pari_printf("K pol: %Ps\n\n", f);
    
    //--------------------------------------------------
    // Define base field K
    struct timespec bnf_started, bnf_finished;
    clock_gettime(CLOCK_MONOTONIC, &bnf_started);
    K = Buchall(f, nf_FORCE, DEFAULTPREC);
    clock_gettime(CLOCK_MONOTONIC, &bnf_finished);
    //K = Buchall_param(f, 1.5, 1.5, 4, nf_FORCE, DEFAULTPREC);
    //--------------------------------------------------
    // Discriminant
    D = nf_get_disc(bnf_get_nf(K));
    pari_printf("Discriminant: %Ps\n\n", D);
    pari_printf("Root discriminant: %Ps\n\n", gsqrtn(gabs(D, DEFAULTPREC), stoi(nf_get_degree(bnf_get_nf(K))), NULL, DEFAULTPREC));
    
    //--------------------------------------------------
    // Check Galois
    if (MY_DEBUGLEVEL >= 1){ my_check_galois(K); }
    
    //--------------------------------------------------        

    // Factor discriminant
    D_prime_vect = gel(factor(D), 1);
    
    //--------------------------------------------------
    // Class group of K (cycle type)
    Kcyc = bnf_get_cyc(K);
    pari_printf("K cyc: %Ps\n\n", Kcyc);
    // pari_printf("r_2(K): %ld\n\n", nf_get_r2(bnf_get_nf(K)));
    //--------------------------------------------------
    // Test if p divides the class number. If not, then H^1(X, Z/pZ) = 0 and there is nothing to compute. 
    my_test_p_rank(K, p_int);
    
    //--------------------------------------------------
    // Find data of unramified extensions
    // my_unramified_p_extensions(K, p, D_prime_vect);
    
    //-------------------------------------------------------------------------------------------------------------------
    // Define polynomials for the generating fields for the part of the Hilbert class field corresp to Cl(K)/p. 
    //-------------------------------------------------------------------------------------------------------------------
    
    //--------------------------------------------------
    // Find generators for the p-torsion of the class group
    struct timespec ja_started, ja_finished;
    clock_gettime(CLOCK_MONOTONIC, &ja_started);
    J_vect = my_find_p_gens(K, p);
    p_rk = lg(J_vect)-1;
    pari_printf("p-rank: %d --> This is the rank of H^1(X,Z/pZ) and H^2(X_fl, mu_p)\n\n", p_rk);
    //--------------------------------------------------

    const char *arithmetic_audit = getenv("MASSEY_ARITHMETIC_AUDIT");
    if (arithmetic_audit && strcmp(arithmetic_audit, "1") == 0)
    {
        GEN audit_units_mod_p = my_find_units_mod_p(K, p);
        Ja_vect =
            my_find_Ja_vect(K, J_vect, p, audit_units_mod_p);
        my_run_arithmetic_audit(K, p, D, Ja_vect, D_prime_vect);
        pari_close();
        return 0;
    }

    const char *rank3_test = getenv("MASSEY_RANK3_TEST");
    if (rank3_test && strcmp(rank3_test, "1") == 0)
    {
        GEN rank3_units_mod_p = my_find_units_mod_p(K, p);
        Ja_vect =
            my_find_Ja_vect(K, J_vect, p, rank3_units_mod_p);
        clock_gettime(CLOCK_MONOTONIC, &ja_finished);
        my_run_rank3_test(
            K, p, D, Kcyc, Ja_vect, D_prime_vect,
            my_elapsed_seconds(bnf_started, bnf_finished),
            my_elapsed_seconds(ja_started, ja_finished),
            start_time);
        pari_close();
        return 0;
    }

    // GEN subgroups = subgrouplist0(bnf_get_cyc(K), stoi(657), 0);

    // Here we generate all subgroups of index p in the class group
    
    GEN subgroups = subgrouplist0(bnf_get_cyc(K), mkvec(p), 0);
    // DEBUG_PRINT(0, "subgroups: %Ps\n\n", subgroups);
    // // pari_close();
    // // exit(0);

    // Here we pick out the subgroups of index p corresponding those p-extensions with smallest class group, but still forming a basis for H^1(X, Z/pZ) 
    GEN best_subgroups = my_best_subgroups(K, p_rk, subgroups, D_prime_vect);
    p_ClFld_pol = bnrclassfield(K, best_subgroups, 0, DEFAULTPREC);

    DEBUG_PRINT(1, "best_subgroups: %Ps\n\n", best_subgroups);
    // p_ClFld_pol = bnrclassfield(K, mkvec2(gel(subgroups, 1),gel(subgroups, 2)), 0, DEFAULTPREC);
    // // p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);
    // DEBUG_PRINT(0, "p Cl Fld (allowing ramification at infinity): %Ps\n\n", p_ClFld_pol);

    // If we don't care which subgroups we use, we can use the default:
    // p_ClFld_pol = bnrclassfield(K, p, 0, DEFAULTPREC);

    DEBUG_PRINT(1, "p Cl Fld: %Ps\n", p_ClFld_pol);
    DEBUG_PRINT(1, ANSI_COLOR_GREEN "Found!\n\n" ANSI_COLOR_RESET);
    
    

    //--------------------------------------------------
    // find generators for the group of units modulo p
    GEN units_mod_p = my_find_units_mod_p(K, p);
    DEBUG_PRINT(1, "Nr of units mod p: %ld\n", glength(units_mod_p));
    //--------------------------------------------------
    // Define r_rk -- the rank of H^2(X, Z/pZ)
    r_rk = glength(J_vect)+glength(units_mod_p);
    pari_printf("r-rank: %d --> This is the rank of H^2(X,Z/pZ) and H^1(X_fl, mu_p)\n\n", r_rk);
    //--------------------------------------------------

    //--------------------------------------------------
    // Define the extensions generating the p-part of the Hilbert class field corresponding to CL(K)/p
    GEN K_ext = my_ext(K, p_ClFld_pol, p, p_rk, D_prime_vect);
    // pari_printf("Extensions found\n\n");
    //--------------------------------------------------

    //--------------------------------------------------
    // Find generators for H^1(X, mu_p), which is dual to H^2(X, Z/pZ)
    Ja_vect = my_find_Ja_vect(K, J_vect, p, units_mod_p);
    //pari_printf("Ja_vect: %Ps\n\n", Ja_vect);
    //--------------------------------------------------

    const char *diagnostics = getenv("MASSEY_DIAGNOSTICS");
    if (diagnostics && strcmp(diagnostics, "1") == 0
        && p_int == 5 && equalii(D, stoi(-90868)))
    {
        my_validate_p_rank_helper();
        long class_p_rank = my_p_class_rank(K, p);
        long discriminant_prime_count = glength(D_prime_vect);
        if (class_p_rank != 2)
            pari_err(e_MISC, "class-group p-rank regression failed");
        pari_printf(
            "  class-group p-rank = %ld; "
            "D_prime_vect length = %ld (independent quantities)\n",
            class_p_rank, discriminant_prime_count);

        GEN family =
            my_secondary_norm_basis_family(
                K, p, Ja_vect, D_prime_vect);
        if (itos(gel(family, 1)) != class_p_rank)
            pari_err(e_MISC,
                     "basis-family dimension is not the class-group p-rank");
        GEN D_a = gmael(family, 2, 1);
        GEN D_b = gmael(family, 2, 2);
        GEN D_ab = gmael(family, 3, 1);
        GEN expected_a =
            mkmat2(mkcol2(stoi(0), stoi(3)),
                   mkcol2(stoi(0), stoi(1)));
        GEN expected_b =
            mkmat2(mkcol2(stoi(0), stoi(0)),
                   mkcol2(stoi(0), stoi(0)));
        GEN expected_ab =
            mkmat2(mkcol2(stoi(2), stoi(3)),
                   mkcol2(stoi(4), stoi(1)));
        if (!gequal(D_a, expected_a)
            || !gequal(D_b, expected_b)
            || !gequal(D_ab, expected_ab))
            pari_err(e_MISC,
                     "secondary norm basis-family regression failed");

        GEN t_c1 = mkcol2(stoi(1), stoi(4));
        GEN t_c2 = mkcol2(stoi(1), stoi(3));
        GEN t_2c1 = mkcol2(stoi(2), stoi(3));

        GEN D_c1 =
            my_secondary_norm_operator(
                K, p, t_c1, Ja_vect, D_prime_vect);
        GEN D_c2 =
            my_secondary_norm_operator(
                K, p, t_c2, Ja_vect, D_prime_vect);
        GEN expected_c1 =
            mkmat2(mkcol2(stoi(3), stoi(3)),
                   mkcol2(stoi(1), stoi(1)));
        GEN expected_c2 =
            mkmat2(mkcol2(stoi(1), stoi(3)),
                   mkcol2(stoi(2), stoi(1)));
        if (!gequal(D_c1, expected_c1))
            pari_err(e_MISC,
                     "prescribed-character regression failed for t=[1,4]");
        if (!gequal(D_c2, expected_c2))
            pari_err(e_MISC,
                     "prescribed-character regression failed for t=[1,3]");

        GEN D_2c1 =
            my_secondary_norm_operator(
                K, p, t_2c1, Ja_vect, D_prime_vect);
        GEN expected_2c1 = FpM_Fp_mul(D_c1, stoi(4), p);
        if (!gequal(D_2c1, expected_2c1))
            pari_err(e_MISC,
                     "prescribed-character scaling regression failed");

        GEN reconstructed_c1 =
            my_reconstruct_secondary_norm(p, family, t_c1);
        GEN reconstructed_c2 =
            my_reconstruct_secondary_norm(p, family, t_c2);
        GEN reconstructed_2c1 =
            my_reconstruct_secondary_norm(p, family, t_2c1);
        if (!gequal(reconstructed_c1, D_c1)
            || !gequal(reconstructed_c2, D_c2)
            || !gequal(reconstructed_2c1, D_2c1))
            pari_err(e_MISC,
                     "quadratic reconstruction regression failed");

        GEN delta_ab =
            my_secondary_norm_delta_basis(p, family, 1, 2);
        GEN delta_ba =
            my_secondary_norm_delta_basis(p, family, 2, 1);
        if (!gequal(delta_ab, delta_ba))
            pari_err(e_MISC, "DeltaD symmetry regression failed");

        GEN T = my_triple_massey_word_matrix(p, family);
        GEN expected_T =
            gp_read_str(
                "[0,3,4,0,3,0,0,0;"
                "0,1,3,0,1,0,0,0]");
        if (!gequal(T, expected_T))
            pari_err(e_MISC,
                     "triple Massey word-matrix regression failed");
        if (FpM_rank(T, p) != 1)
            pari_err(e_MISC,
                     "triple Massey word-matrix rank is not 1");
        my_validate_triple_massey_identities(p, T, 2);

        GEN characters[2] = { t_c1, t_c2 };
        GEN operators[2] = { D_c1, D_c2 };
        for (long x_index = 0; x_index < 2; ++x_index)
            for (long y_index = 0; y_index < 2; ++y_index)
            {
                GEN Q =
                    my_evaluate_character_on_operator(
                        p, characters[y_index],
                        operators[x_index]);
                GEN contracted =
                    my_triple_massey_contract(
                        p, T, 2,
                        characters[x_index],
                        characters[x_index],
                        characters[y_index]);
                if (!gequal(Q, contracted))
                    pari_err(e_MISC,
                             "Q(x,y) tensor contraction regression failed");
                pari_printf(
                    "  Q check x=%Ps, y=%Ps: direct=%Ps, tensor=%Ps, OK\n",
                    gtovec(characters[x_index]),
                    gtovec(characters[y_index]),
                    gtovec(Q), gtovec(contracted));
            }

        pari_printf(
            "\nMASSEY_DIAGNOSTICS quadratic/tensor regressions OK\n");
        pari_printf("  D_a = %Ps\n", D_a);
        pari_printf("  D_b = %Ps\n", D_b);
        pari_printf("  D_[a+b] = %Ps\n", D_ab);
        pari_printf("  D_[1,4] = %Ps\n", D_c1);
        pari_printf("  D_[1,3] = %Ps\n", D_c2);
        pari_printf("  D_[2,3] = %Ps = 4 * D_[1,4]\n", D_2c1);
        pari_printf("  reconstructed D_[1,4] = %Ps\n", reconstructed_c1);
        pari_printf("  reconstructed D_[1,3] = %Ps\n", reconstructed_c2);
        pari_printf("  reconstructed D_[2,3] = %Ps\n", reconstructed_2c1);
        pari_printf("  DeltaD(a,b) = %Ps\n", delta_ab);
        pari_printf("  triple Massey evaluation word matrix T = %Ps\n", T);
        pari_printf("  rank(T) = %ld\n", FpM_rank(T, p));
        pari_printf("  shuffle identities = OK\n");
    }
    
    //--------------------------------------------------
    // CUP PRODUCTS
    //--------------------------------------------------
    // Defines a matrix over F_p with index (i*k, j) corresponding to 
    // < x_i\cup x_k, (a_j, J_j) > if i is not equal to j and
    // < B(x_i), (a_j, J_j) > if i=j. 
    // Here < - , - > denotes the Artin--Verdier pairing, which may be computed using our cup product formula and the Artin symbol. 
    //--------------------------------------------------

    //--------------------------------------------------
    // Non-parallel version
    // int mat_rk = my_relations(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk);
    //--------------------------------------------------
    // Parallell computation of the cup products. These are always zero for imaginary quadratic fields. 
    // int mat_rk = my_relations_par(K_ext, K, p, p_rk, Ja_vect, r_rk);
    //---------------------

    //--------------------------------------------------
    // HIGHER MASSEY PRODUCTS (This is only implemented for some Massey products of the form < x, x, ..., x, y > as seen below, but still very useful)
    //--------------------------------------------------
    // Defines a matrix over F_p with index (i*k, j) corresponding to 
    // < x_i, x_i, ..., x_i, x_k, (a_j, J_j) > if i is not equal to j and
    //--------------------------------------------------
    int mat_rk = 0;
    if ((mat_rk<3 && p_int>2) || (mat_rk==0 && p_int==2))
    {
        my_print_massey(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk, best_subgroups);
    }    
    
    //--------------------------------------------------

    DEBUG_PRINT(0, ANSI_COLOR_GREEN "Done! \n \n" ANSI_COLOR_YELLOW);
   
    pari_close();






    //--------------------------------------------------
    // Compute the CPU time
    clock_t duration = (clock()-start) / 1000; // Compute CPU duration in microseconds

    // Compute the actual time
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    
    // Compute actual duration in milliseconds
    long duration_ns = (end_time.tv_sec - start_time.tv_sec) * 1e9 + (end_time.tv_nsec - start_time.tv_nsec);
    long duration_ms = duration_ns / 1e6; // Convert to milliseconds

    // Convert to minutes, seconds, and milliseconds
    msec = duration_ms % 1000;
    sec = (duration_ms / 1000) % 60;
    min = duration_ms / 60000;

    // Print actual elapsed time
    printf(ANSI_COLOR_YELLOW "Actual time: %d min, %d sec, %d msec\n\n" ANSI_COLOR_RESET, min, sec, msec);
    
    // Compute the CPU time
    msec = duration%1000000;
    sec = (duration/1000)%60;
    min = duration/60000;

    printf (ANSI_COLOR_YELLOW "CPU time: %d min, %d,%d sec\n" ANSI_COLOR_RESET, min, sec, msec);
    //--------------------------------------------------
    printf("\n---------------------------------------------------------------------------------------------------------\nEnd program\n---------------------------------------------------------------------------------------------------------\n\n");
    return 0;
}
