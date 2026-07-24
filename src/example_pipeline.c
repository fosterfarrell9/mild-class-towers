// MIT License

/**
 * @file example_pipeline.c
 * @brief Audited orchestration for compact rank-three example records.
 *
 * Arithmetic candidates come from the existing relative-BNF implementation
 * and are accepted only after its exact audit. Tensor reconstruction and
 * strong-freeness analysis then operate solely on finite-field matrices.
 */

#include <stdio.h>
#include <pari/pari.h>
#include "../headers/example_pipeline.h"
#include "../headers/massey_tensor.h"
#include "../headers/relation_degree3.h"
#include "../headers/secondary_norm.h"

static GEN
example_field(const char *name, GEN value)
{
    return mkvec2(strtoGENstr(name), value);
}

static void
example_write_record(const char *path, GEN record)
{
    FILE *file = fopen(path, "wb");
    if (!file) pari_err_FILE("example result", path);
    pari_fprintf(file, "%Ps\n", record);
    if (fclose(file) != 0) pari_err_FILE("example result", path);
}

GEN
my_compute_example_result(
    GEN K, GEN p, GEN polynomial, GEN discriminant,
    GEN Ja_vect, GEN D_prime_vect, const char *output_path,
    long strong_search_limit)
{
    pari_sp av = avma;
    if (my_p_class_rank(K, p) != 3)
        pari_err(e_MISC, "example pipeline requires p-class rank 3");
    if (glength(Ja_vect) != 3)
        pari_err(e_MISC, "example pipeline requires three arithmetic inputs");
    if (bnfcertify0(K, 0) != 1)
        pari_err(e_MISC, "example pipeline base BNF certification failed");

    pari_printf("BASE_BNF_CERTIFIED\n");
    pari_printf("P_CLASS_RANK 3\n");

    /* Require exact acceptance checks for every computed matrix column. */
    my_secondary_norm_require_exact_audit(1);
    GEN family =
        my_secondary_norm_basis_family(K, p, Ja_vect, D_prime_vect);
    my_secondary_norm_require_exact_audit(0);
    pari_printf("SECONDARY_NORMS VERIFIED\n");

    GEN basis = gel(family, 2);
    GEN pairs = gel(family, 3);
    GEN samples = cgetg(7, t_VEC);
    for (long i = 1; i <= 3; ++i) gel(samples, i) = gel(basis, i);
    for (long i = 1; i <= 3; ++i) gel(samples, i + 3) = gel(pairs, i);

    GEN doubled = cgetg(4, t_VEC);
    for (long i = 1; i <= 3; ++i)
    {
        GEN t = zerocol(3);
        gel(t, i) = stoi(2);
        GEN reconstructed = my_reconstruct_secondary_norm(p, family, t);
        GEN expected = FpM_Fp_mul(gel(basis, i), stoi(4), p);
        if (!gequal(reconstructed, expected))
            pari_err(e_MISC, "example pipeline quadratic scaling failed");
        gel(doubled, i) = reconstructed;
    }

    /*
     * Columns are X_i X_j X_k with k varying fastest, then j, then i.
     * The identity check covers outer symmetry, cyclic shuffle, and diagonal
     * vanishing before the matrix is passed to relation analysis.
     */
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
        (q * q * q - 1)
        * (q * q * q - q)
        * (q * q * q - q * q);
    long effective_limit =
        strong_search_limit < 0 ? gl_order : strong_search_limit;
    GEN fields = cgetg(21, t_VEC);
    gel(fields, 1) = example_field("format_version", gen_1);
    gel(fields, 2) = example_field("status", strtoGENstr(status));
    gel(fields, 3) = example_field("p", p);
    gel(fields, 4) = example_field("base_polynomial", polynomial);
    gel(fields, 5) = example_field("base_discriminant", discriminant);
    gel(fields, 6) = example_field("class_group_invariants", bnf_get_cyc(K));
    gel(fields, 7) = example_field("class_number", bnf_get_no(K));
    gel(fields, 8) = example_field("class_group_generators", bnf_get_gen(K));
    gel(fields, 9) = example_field("character_basis_columns", matid(3));
    gel(fields, 10) = example_field("secondary_norm_samples", samples);
    gel(fields, 11) = example_field("doubled_character_checks", doubled);
    gel(fields, 12) = example_field("arithmetic_exact_audit", gen_1);
    gel(fields, 13) = example_field(
        "word_order", strtoGENstr("X_i X_j X_k; k fastest"));
    gel(fields, 14) = example_field("cubic_relation_matrix", T);
    gel(fields, 15) = example_field("cubic_rank", stoi(rank));
    gel(fields, 16) = example_field("strong_freeness_witness", witness);
    gel(fields, 17) = example_field(
        "strong_freeness_candidate_limit", stoi(effective_limit));
    gel(fields, 18) = example_field(
        "strong_freeness_exhaustive_gl",
        strong_search_limit < 0 ? gen_1 : gen_0);
    gel(fields, 19) = example_field(
        "MILD", strtoGENstr(mild ? "PROVED" : "UNKNOWN"));
    gel(fields, 20) = example_field(
        "CD", mild ? gen_2 : strtoGENstr("UNKNOWN"));

    GEN record = gerepilecopy(av, fields);
    example_write_record(output_path, record);
    return record;
}
