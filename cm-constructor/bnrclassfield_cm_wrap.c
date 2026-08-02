/* Link-time experiment wrapper: replace only the class-field search oracle
 * with an already recognized CM quintic.  The caller, all exact audits, and
 * every downstream routine remain the versioned pipeline objects.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>

GEN __real_bnrclassfield(GEN bnr, GEN subgroup, long flag, long prec);

static void
wrapper_error(const char *message)
{
    pari_err(e_MISC, "CM bnrclassfield wrapper: %s", message);
}

static GEN
read_expression(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (!file) wrapper_error("cannot open MASSEY_CM_OUTPUT");
    if (fseek(file, 0, SEEK_END) != 0)
        wrapper_error("cannot seek MASSEY_CM_OUTPUT");
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0)
        wrapper_error("cannot size MASSEY_CM_OUTPUT");
    char *text = malloc((size_t)size + 1);
    if (!text) wrapper_error("cannot allocate input buffer");
    if (fread(text, 1, (size_t)size, file) != (size_t)size)
        wrapper_error("cannot read MASSEY_CM_OUTPUT");
    fclose(file);
    text[size] = '\0';
    GEN value = gp_read_str(text);
    free(text);
    return value;
}

GEN
__wrap_bnrclassfield(GEN bnr, GEN subgroup, long flag, long prec)
{
    const char *path = getenv("MASSEY_CM_OUTPUT");
    const char *index_text = getenv("MASSEY_CM_CHARACTER");
    if (!path || !*path)
        return __real_bnrclassfield(bnr, subgroup, flag, prec);
    if (!index_text || !*index_text)
        wrapper_error("MASSEY_CM_CHARACTER is not set");

    long index = atol(index_text);
    if (index < 1 || index > 6)
        wrapper_error("MASSEY_CM_CHARACTER is outside 1..6");

    GEN cm = read_expression(path);
    if (typ(cm) != t_VEC || lg(cm) != 12 || !equaliu(gel(cm, 1), 1))
        wrapper_error("invalid CM output schema");
    if (!gequal(nf_get_pol(bnf_get_nf(bnr)), gel(cm, 4))
        || !gequal(nf_get_disc(bnf_get_nf(bnr)), gel(cm, 5)))
        wrapper_error("CM output belongs to a different base field");

    static const long expected[6][3] = {
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
    };
    GEN entries = gel(cm, 9);
    for (long i = 1; i < lg(entries); ++i)
    {
        GEN entry = gel(entries, i);
        GEN character = gel(entry, 2);
        int match = 1;
        for (long coordinate = 1; coordinate <= 3; ++coordinate)
            if (!equaliu(
                    gel(character, coordinate),
                    expected[index - 1][coordinate - 1]))
            {
                match = 0;
                break;
            }
        if (!match) continue;
        GEN polynomial = gel(entry, 4);
        pari_printf(
            "CM_ORACLE_INJECTED character=%ld polynomial=%Ps "
            "bnrclassfield_Kzeta5=ELIMINATED\n",
            index, polynomial);
        return gcopy(polynomial);
    }
    wrapper_error("requested character is absent from CM output");
    return NULL;
}
