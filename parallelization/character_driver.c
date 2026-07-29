// MIT License

/**
 * @file character_driver.c
 * @brief Compute one prescribed-character secondary norm in its own process.
 *
 * This driver reuses the existing, unmodified modules of Massey-pari: it
 * replicates the deterministic input preparation of the audited pipeline
 * (base BNF, certification, torsion representatives, unit generators,
 * Ahlqvist--Carlson pairs) and then computes the secondary-norm matrix of a
 * single prescribed character with the exact audit enabled.  Running six
 * such processes, one per character, parallelizes a field's arithmetic
 * without touching the proven sequential code path.
 *
 * Usage:
 *   character_driver <p> <polynomial> <index 1..6> <matrix-output>
 *
 * With MASSEY_CERTIFICATE_EXPORT set, the process writes a well-formed
 * partial certificate containing this character's three entries.
 */

#include <stdio.h>
#include <stdlib.h>
#include <pari/pari.h>
#include "../headers/misc_functions.h"
#include "../headers/secondary_norm.h"
#include "../headers/find_cup_matrix.h"

static const long CHARACTER_COORDS[6][3] = {
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
    {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
};

int
main(int argc, char *argv[])
{
    if (argc != 5)
    {
        fprintf(stderr,
                "usage: %s <p> <polynomial> <index 1..6> <matrix-output>\n",
                argv[0]);
        return 2;
    }

    /* Initialization as in the sequential main, with an optionally
     * reduced stack ceiling so several processes fit into the WSL
     * memory cap (MASSEY_PARISTACK_MAX in bytes, default 8 GiB). */
    setvbuf(stdout, NULL, _IOLBF, 0);
    long long stack_max = 1LL << 33;
    const char *cap = getenv("MASSEY_PARISTACK_MAX");
    if (cap && *cap) stack_max = atoll(cap);
    entree ep = {"_worker", 0, (void *)compute_my_relations, 20, "LG", ""};
    pari_init_opts(1L << 30, 1048576,
                   INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    pari_add_function(&ep);
    pari_mt_init();
    paristack_setsize(1L << 30, (size_t)stack_max);
    sd_threadsizemax("2147483648", 0);
    setalldebug(0);

    GEN p = gp_read_str(argv[1]);
    GEN f = gp_read_str(argv[2]);
    long index = atol(argv[3]);
    if (index < 1 || index > 6)
        pari_err(e_MISC, "character index must be in 1..6");

    GEN K = Buchall(f, nf_FORCE, DEFAULTPREC);
    GEN D = nf_get_disc(bnf_get_nf(K));
    pari_printf("CHARACTER_DRIVER index=%ld D=%Ps\n", index, D);

    if (nf_get_degree(bnf_get_nf(K)) != 2
        || nf_get_r2(bnf_get_nf(K)) != 1)
        pari_err(e_MISC, "driver requires an imaginary quadratic field");
    if (my_p_class_rank(K, p) != 3)
        pari_err(e_MISC, "driver requires p-class rank 3");
    if (bnfcertify0(K, 0) != 1)
        pari_err(e_MISC, "base BNF certification failed");
    pari_printf("BASE_BNF_CERTIFIED\n");

    /* Deterministic input preparation, exactly as in the sequential main. */
    GEN D_prime_vect = gel(factor(D), 1);
    GEN J_vect = my_find_p_gens(K, p);
    GEN units_mod_p = my_find_units_mod_p(K, p);
    GEN Ja_vect = my_find_Ja_vect(K, J_vect, p, units_mod_p);
    if (glength(Ja_vect) != 3)
        pari_err(e_MISC, "driver requires three arithmetic inputs");

    GEN t = mkcol3(
        stoi(CHARACTER_COORDS[index - 1][0]),
        stoi(CHARACTER_COORDS[index - 1][1]),
        stoi(CHARACTER_COORDS[index - 1][2]));
    pari_printf("prescribed t = %Ps\n", gtovec(t));

    my_secondary_norm_require_exact_audit(1);
    GEN Dt = my_secondary_norm_operator(K, p, t, Ja_vect, D_prime_vect);
    my_secondary_norm_require_exact_audit(0);

    FILE *out = fopen(argv[4], "wb");
    if (!out) pari_err_FILE("matrix output", argv[4]);
    pari_fprintf(out, "[%Ps,%Ps]\n", t, Dt);
    if (fclose(out) != 0) pari_err_FILE("matrix output", argv[4]);

    pari_printf("CHARACTER_DONE index=%ld D_t=%Ps\n", index, Dt);
    pari_close();
    return 0;
}
