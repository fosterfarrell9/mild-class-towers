// MIT License

/**
 * @file finish_driver.c
 * @brief Assemble a result record from six per-character matrices.
 *
 * Reads the six matrices produced by character_driver, rebuilds the
 * quadratic family, and performs exactly the finite-field part of the
 * sequential pipeline (quadratic-scaling checks, tensor identities, rank,
 * strong-freeness witness search) using the existing, unmodified modules.
 * The written record has the same schema as my_compute_example_result.
 *
 * Usage:
 *   finish_driver <p> <polynomial> <result-output> <limit|exhaustive> \
 *                 <matrix-file-1> ... <matrix-file-6>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pari/pari.h>
#include "../headers/misc_functions.h"
#include "../headers/secondary_norm.h"
#include "../headers/massey_tensor.h"
#include "../headers/relation_degree3.h"
#include "../headers/example_pipeline.h"
#include "../headers/find_cup_matrix.h"

static const long CHARACTER_COORDS[6][3] = {
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
    {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
};

static GEN
record_field(const char *name, GEN value)
{
    return mkvec2(strtoGENstr(name), value);
}

int
main(int argc, char *argv[])
{
    if (argc != 11)
    {
        fprintf(stderr,
                "usage: %s <p> <polynomial> <result-output> "
                "<limit|exhaustive> <matrix-file-1> ... <matrix-file-6>\n",
                argv[0]);
        return 2;
    }

    entree ep = {"_worker", 0, (void *)compute_my_relations, 20, "LG", ""};
    pari_init_opts(1L << 30, 1048576,
                   INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    pari_add_function(&ep);
    pari_mt_init();
    paristack_setsize(1L << 30, 1L << 33);
    sd_threadsizemax("2147483648", 0);
    setalldebug(0);

    GEN p = gp_read_str(argv[1]);
    GEN f = gp_read_str(argv[2]);
    const char *output_path = argv[3];
    long strong_search_limit;
    if (strcmp(argv[4], "exhaustive") == 0)
        strong_search_limit = -1;
    else
    {
        char *end = NULL;
        strong_search_limit = strtol(argv[4], &end, 10);
        if (!argv[4][0] || !end || *end || strong_search_limit <= 0)
            pari_err(e_MISC, "invalid strong-search limit");
    }

    GEN K = Buchall(f, nf_FORCE, DEFAULTPREC);
    GEN D = nf_get_disc(bnf_get_nf(K));
    if (my_p_class_rank(K, p) != 3)
        pari_err(e_MISC, "finish driver requires p-class rank 3");
    if (bnfcertify0(K, 0) != 1)
        pari_err(e_MISC, "finish driver base BNF certification failed");
    pari_printf("BASE_BNF_CERTIFIED\n");
    pari_printf("P_CLASS_RANK 3\n");

    /* Collect the six matrices; verify the stored characters match. */
    GEN basis = cgetg(4, t_VEC);
    GEN pairs = cgetg(4, t_VEC);
    for (long q = 0; q < 6; ++q)
    {
        GEN pair = gp_read_file(argv[5 + q]);
        if (typ(pair) != t_VEC || glength(pair) != 2)
            pari_err(e_MISC, "matrix file %s is not [t, D_t]", argv[5 + q]);
        GEN t = gel(pair, 1);
        GEN expected = mkcol3(
            stoi(CHARACTER_COORDS[q][0]),
            stoi(CHARACTER_COORDS[q][1]),
            stoi(CHARACTER_COORDS[q][2]));
        if (!gequal(t, expected))
            pari_err(e_MISC, "matrix file %s has the wrong character",
                     argv[5 + q]);
        if (q < 3) gel(basis, q + 1) = gel(pair, 2);
        else gel(pairs, q - 2) = gel(pair, 2);
    }
    GEN family = mkvec3(stoi(3), basis, pairs);
    pari_printf("SECONDARY_NORMS COLLECTED\n");

    GEN samples = cgetg(7, t_VEC);
    for (long i = 1; i <= 3; ++i) gel(samples, i) = gel(basis, i);
    for (long i = 1; i <= 3; ++i) gel(samples, i + 3) = gel(pairs, i);

    /* Finite algebra identical to my_compute_example_result. */
    GEN doubled = cgetg(4, t_VEC);
    for (long i = 1; i <= 3; ++i)
    {
        GEN t = zerocol(3);
        gel(t, i) = stoi(2);
        GEN reconstructed = my_reconstruct_secondary_norm(p, family, t);
        GEN expected = FpM_Fp_mul(gel(basis, i), stoi(4), p);
        if (!gequal(reconstructed, expected))
            pari_err(e_MISC, "quadratic scaling failed");
        gel(doubled, i) = reconstructed;
    }

    GEN T = my_triple_massey_word_matrix(p, family);
    my_validate_triple_massey_identities(p, T, 3);
    long rank = FpM_rank(T, p);
    pari_printf("MASSEY_RANK %ld\n", rank);

    GEN witness = rank == 3
        ? my_find_strongly_free_witness(T, p, strong_search_limit)
        : cgetg(1, t_VEC);
    int mild = rank == 3 && glength(witness) != 0;
    const char *status = rank < 3
        ? MASSEY_EXAMPLE_STATUS_RANK_LT_3
        : mild ? MASSEY_EXAMPLE_STATUS_PROVED
               : MASSEY_EXAMPLE_STATUS_NO_WITNESS;
    if (mild)
    {
        pari_printf("STRONG_FREENESS PASS\n");
        pari_printf("LEADING_WORDS %Ps\n", gel(witness, 4));
        pari_printf("MILD PASS\n");
    }
    else
        pari_printf("MILD UNKNOWN\n");

    long q = itos(p);
    long gl_order =
        (q * q * q - 1) * (q * q * q - q) * (q * q * q - q * q);
    long effective_limit =
        strong_search_limit < 0 ? gl_order : strong_search_limit;
    GEN fields = cgetg(21, t_VEC);
    gel(fields, 1) = record_field("format_version", gen_1);
    gel(fields, 2) = record_field("status", strtoGENstr(status));
    gel(fields, 3) = record_field("p", p);
    gel(fields, 4) = record_field("base_polynomial", f);
    gel(fields, 5) = record_field("base_discriminant", D);
    gel(fields, 6) = record_field("class_group_invariants", bnf_get_cyc(K));
    gel(fields, 7) = record_field("class_number", bnf_get_no(K));
    gel(fields, 8) = record_field("class_group_generators", bnf_get_gen(K));
    gel(fields, 9) = record_field("character_basis_columns", matid(3));
    gel(fields, 10) = record_field("secondary_norm_samples", samples);
    gel(fields, 11) = record_field("doubled_character_checks", doubled);
    gel(fields, 12) = record_field("arithmetic_exact_audit", gen_1);
    gel(fields, 13) = record_field(
        "word_order", strtoGENstr("X_i X_j X_k; k fastest"));
    gel(fields, 14) = record_field("cubic_relation_matrix", T);
    gel(fields, 15) = record_field("cubic_rank", stoi(rank));
    gel(fields, 16) = record_field("strong_freeness_witness", witness);
    gel(fields, 17) = record_field(
        "strong_freeness_candidate_limit", stoi(effective_limit));
    gel(fields, 18) = record_field(
        "strong_freeness_exhaustive_gl",
        strong_search_limit < 0 ? gen_1 : gen_0);
    gel(fields, 19) = record_field(
        "MILD", strtoGENstr(mild ? "PROVED" : "UNKNOWN"));
    gel(fields, 20) = record_field(
        "CD", mild ? gen_2 : strtoGENstr("UNKNOWN"));

    FILE *file = fopen(output_path, "wb");
    if (!file) pari_err_FILE("result output", output_path);
    pari_fprintf(file, "%Ps\n", fields);
    if (fclose(file) != 0) pari_err_FILE("result output", output_path);

    pari_printf("RESULT_WRITTEN %s\n", output_path);
    pari_close();
    return 0;
}
