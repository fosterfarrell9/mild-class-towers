// MIT License

/**
 * @file upgrade_certificate.c
 * @brief Rewrite a format 1 certificate as format 2 by recording its bases.
 *
 * A format 1 certificate stores elements and ideals as coordinates with
 * respect to an integral basis it does not name.  Since PARI's integral basis
 * is LLL-reduced and not canonical, such a certificate is only meaningful on
 * a machine whose nfinit reproduces the basis its generator used.
 *
 * This program supplies what is missing.  It must run on a machine where the
 * certificate verifies -- there and only there does nfinit return the basis
 * the coordinates were written against -- and it recomputes nothing else: the
 * arithmetic of the certificate is copied through unchanged.
 *
 * Usage:
 *   upgrade_certificate <format-1 certificate> <output>
 *
 * Afterwards check the result with the verifier, on this machine and on
 * another one; only the second run tests what format 2 is for.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>

enum {
    CERT_FORMAT = 1,
    CERT_PARI_VERSION,
    CERT_P,
    CERT_BASE_POLYNOMIAL,
    CERT_DISCRIMINANT,
    CERT_BASE_DATA,
    CERT_ENTRIES
};

enum { ENTRY_ABSOLUTE_POLYNOMIAL = 5 };

static void
die(const char *message)
{
    pari_fprintf(stderr, "upgrade_certificate: %s\n", message);
    pari_close();
    exit(EXIT_FAILURE);
}

static GEN
read_file(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (!file) die("cannot open certificate");
    if (fseek(file, 0, SEEK_END) != 0) die("cannot seek certificate");
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0)
        die("cannot size certificate");
    char *text = malloc((size_t)size + 1);
    if (!text) die("cannot allocate certificate buffer");
    if (fread(text, 1, (size_t)size, file) != (size_t)size)
        die("cannot read certificate");
    fclose(file);
    text[size] = '\0';
    GEN parsed = gp_read_str(text);
    free(text);
    return parsed;
}

int
main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "usage: %s <format-1 certificate> <output>\n",
                argv[0]);
        return 2;
    }

    pari_init_opts(1L << 30, 1048576,
                   INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    paristack_setsize(1L << 30, 1L << 33);

    GEN certificate = read_file(argv[1]);
    if (typ(certificate) != t_VEC || lg(certificate) != 8)
        die("invalid top-level certificate schema");
    if (!equaliu(gel(certificate, CERT_FORMAT), 1))
        die("input is not a format 1 certificate");
    if (!equaliu(gel(certificate, CERT_PARI_VERSION), PARI_VERSION_CODE))
        die("PARI version differs from the certificate generator");

    /* The base field, for its integral basis.  bnfcertify is the verifier's
     * business, not ours; here the field is only a source of a basis. */
    GEN K = Buchall(gel(certificate, CERT_BASE_POLYNOMIAL),
                    nf_FORCE, DEFAULTPREC);
    if (!gequal(nf_get_disc(bnf_get_nf(K)), gel(certificate, CERT_DISCRIMINANT)))
        die("base discriminant mismatch");

    GEN base_data = gel(certificate, CERT_BASE_DATA);
    if (typ(base_data) != t_VEC || lg(base_data) != 6)
        die("invalid base metadata");

    FILE *out = fopen(argv[2], "wb");
    if (!out) die("cannot open output");

    pari_fprintf(out, "[2,%Ps,%Ps,%Ps,%Ps,[%Ps,%Ps,%Ps,%Ps,%Ps,%Ps],\n[",
                 gel(certificate, CERT_PARI_VERSION), gel(certificate, CERT_P),
                 gel(certificate, CERT_BASE_POLYNOMIAL),
                 gel(certificate, CERT_DISCRIMINANT),
                 gel(base_data, 1), gel(base_data, 2), gel(base_data, 3),
                 gel(base_data, 4), gel(base_data, 5),
                 nf_get_zk(bnf_get_nf(K)));

    GEN entries = gel(certificate, CERT_ENTRIES);
    if (typ(entries) != t_VEC) die("invalid entry list");

    /* The eighteen entries share six class fields, three entries each, so
     * cache by defining polynomial rather than run nfinit eighteen times. */
    GEN seen_polynomials = cgetg(lg(entries), t_VEC);
    GEN seen_bases = cgetg(lg(entries), t_VEC);
    long cached = 0;

    for (long index = 1; index < lg(entries); ++index)
    {
        GEN entry = gel(entries, index);
        if (typ(entry) != t_VEC || lg(entry) != 14)
            die("invalid entry schema");
        GEN polynomial = gel(entry, ENTRY_ABSOLUTE_POLYNOMIAL);

        GEN basis = NULL;
        for (long i = 1; i <= cached; ++i)
            if (gequal(gel(seen_polynomials, i), polynomial))
            {
                basis = gel(seen_bases, i);
                break;
            }
        if (!basis)
        {
            basis = nf_get_zk(nfinit0(polynomial, 0, DEFAULTPREC));
            ++cached;
            gel(seen_polynomials, cached) = polynomial;
            gel(seen_bases, cached) = basis;
        }

        if (index > 1) fputs(",\n", out);
        pari_fprintf(out, "[");
        for (long field = 1; field < lg(entry); ++field)
        {
            /* %Ps prints a t_STR bare, which would read back as a variable
             * name rather than the character label. */
            if (typ(gel(entry, field)) == t_STR)
                pari_fprintf(out, "\"%Ps\",", gel(entry, field));
            else
                pari_fprintf(out, "%Ps,", gel(entry, field));
        }
        pari_fprintf(out, "%Ps]", basis);
        fflush(out);
    }

    pari_fprintf(out, "\n]]\n");
    if (fclose(out) != 0) die("cannot close output");
    pari_printf("upgraded to format 2: %s (%ld entries, %ld class fields)\n",
                argv[2], lg(entries) - 1, cached);
    pari_close();
    return 0;
}
