#define _POSIX_C_SOURCE 200809L

/*
 * Deterministic CM constructor for the six cyclic quintic subfields used by
 * Massey-pari.  This is experiment-only code.  The double-eta evaluation
 * below is adapted from PARI/GP 2.17.4 src/modules/stark.c (GPL-2+); in
 * particular it preserves Schertz's eta-multiplier corrections.
 */

#include <pari/pari.h>
#include <pari/paripriv.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const long CHARACTER_COORDS[6][3] = {
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
    {1, 1, 0}, {1, 0, 1}, {0, 1, 1},
};

static const char *CHARACTER_LABELS[6] = {
    "a", "b", "c", "a+b", "a+c", "b+c",
};

struct gpq_data {
    long p, q;
    GEN sqd;
    GEN u, D;
    GEN pq, pq2;
    GEN qfpq;
};

static double
monotonic_seconds(void)
{
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return 0.0;
    return (double)now.tv_sec + (double)now.tv_nsec / 1.0e9;
}

static long
gcd24(long n)
{
    static const long values[24] = {
        24, 1, 2, 3, 4, 1, 6, 1, 8, 3, 2, 1,
        12, 1, 2, 3, 8, 1, 6, 1, 4, 3, 2, 1,
    };
    return values[n % 24];
}

static int
has_exponent_two(GEN form)
{
    GEN a = gel(form, 1), b = gel(form, 2), c = gel(form, 3);
    return !signe(b) || absequalii(a, b) || equalii(a, c);
}

/* PARI's deterministic choice used by quadhilbertimag. */
static void
initialize_double_eta(GEN D, struct gpq_data *data)
{
    const long capacity = 6547;
    const ulong max_q = 50000;
    GEN primes = cgetg(capacity + 1, t_VECSMALL);
    GEN order_two_forms = cgetg(capacity + 1, t_VEC);
    GEN gcds = cgetg(capacity + 1, t_VECSMALL);
    forprime_t iterator;
    long length = 1;
    double best = 0.0;
    ulong q;

    u_forprime_init(&iterator, 2, ULONG_MAX);
    data->D = D;
    data->p = data->q = 0;
    for (;;)
    {
        GEN Q;
        long i, gcd_q, modulus;
        int order_two, store;
        double quotient;

        q = u_forprime_next(&iterator);
        if (best > 0.0 && q >= max_q) break;
        if (kroiu(D, q) < 0) continue;
        Q = qfbred_i(primeform_u(D, q));
        if (is_pm1(gel(Q, 1))) continue;

        store = 1;
        order_two = has_exponent_two(Q);
        gcds[length] = gcd_q = gcd24((long)q - 1);
        modulus = 24 / gcd_q;
        primes[length] = (long)q;
        gel(order_two_forms, length) = order_two ? Q : NULL;
        quotient = (q + 1) / (double)(q - 1);

        for (i = 1; i < length; ++i)
        {
            long p = primes[i], gcd_p = gcds[i];
            double candidate;
            if (order_two && gel(order_two_forms, i)
                && !gequal(gel(order_two_forms, i), Q))
                continue;
            if (gcd_p % gcd_q == 0) store = 0;
            if ((p - 1) % modulus) continue;
            candidate = quotient * (p + 1) / (double)(p - 1);
            if (candidate > best)
            {
                store = 0;
                best = candidate;
                data->q = (long)q;
                data->p = p;
            }
            if (best > 0.0) break;
        }
        if (store && quotient * quotient > best)
            if (++length >= capacity)
                pari_err_BUG("CM constructor: eta prime capacity");
        if (best == 0.0 && gcd_q >= 12 && umodiu(D, q))
        {
            double candidate = quotient * q / (double)(q - 1);
            if (candidate > best)
            {
                best = candidate;
                data->q = data->p = (long)q;
            }
        }
        if ((primes[1] + 1) * quotient <= (primes[1] - 1) * best)
            break;
    }
    if (!data->p || !data->q)
        pari_err_BUG("CM constructor: no double-eta primes");
}

static void
finish_double_eta_initialization(struct gpq_data *data)
{
    GEN prime_p_form = primeform_u(data->D, data->p), u;
    data->pq = muluu((ulong)data->p, (ulong)data->q);
    data->pq2 = shifti(data->pq, 1);
    if (data->p == data->q)
    {
        GEN square = qfbcompraw(prime_p_form, prime_p_form);
        u = gel(square, 2);
        data->u = modii(u, data->pq2);
        data->qfpq = qfbred_i(square);
    }
    else
    {
        GEN prime_q_form = primeform_u(data->D, data->q);
        data->u = Z_chinese(
            gel(prime_p_form, 2), gel(prime_q_form, 2),
            utoipos((ulong)data->p << 1),
            utoipos((ulong)data->q << 1));
        data->qfpq = qfbcomp_i(prime_p_form, prime_q_form);
    }
}

/* Evaluate the correctly corrected double-eta class invariant for every
 * form.  PARI's full-class-polynomial routine skips the conjugate partner;
 * here both partners are needed because traces are grouped by kernel coset.
 */
static GEN
double_eta_value(GEN form, struct gpq_data *data)
{
    pari_sp av = avma;
    long a = itos(gel(form, 1));
    long b = itos(gel(form, 2));
    long c = itos(gel(form, 3));
    long p = data->p, q = data->q;
    GEN w, value;

    if (p == 2 && a % q == 0 && (a & b & 1) && !(c & 1))
    {
        lswap(a, c);
        b = -b;
    }
    if (a % p == 0 || a % q == 0)
    {
        while (c % p == 0 || c % q == 0)
        {
            c += a + b;
            b += a << 1;
        }
        lswap(a, c);
        b = -b;
    }
    w = Z_chinese(
        data->u, stoi(-b), data->pq2, utoipos((ulong)a << 1));
    value = double_eta_quotient(
        utoipos((ulong)a), w, data->D, p, q, data->pq, data->sqd);
    return gerepileupto(av, value);
}

static GEN
form_from_base_ideal(GEN ideal, GEN D)
{
    GEN a, residue, b, c, numerator;
    if (typ(ideal) != t_MAT || lg(ideal) != 3
        || lg(gel(ideal, 1)) != 3 || lg(gel(ideal, 2)) != 3)
        pari_err_TYPE("CM constructor: quadratic ideal HNF", ideal);
    if (!gequal0(gcoeff(ideal, 2, 1))
        || !gequal1(gcoeff(ideal, 2, 2)))
        pari_err_BUG("CM constructor: unexpected quadratic ideal HNF");
    a = gcoeff(ideal, 1, 1);
    residue = gcoeff(ideal, 1, 2);
    b = subii(negi(shifti(residue, 1)), gen_1);
    numerator = subii(sqri(b), D);
    c = diviiexact(numerator, shifti(a, 2));
    return qfbred_i(mkqfb(a, b, c, D));
}

static GEN
principal_form(GEN D)
{
    GEN b = mpodd(D) ? gen_1 : gen_0;
    GEN c = shifti(subii(sqri(b), D), -2);
    return mkqfb(gen_1, b, c, D);
}

static GEN
form_powers(GEN generator, long order, GEN identity)
{
    GEN powers = cgetg(order + 1, t_VEC);
    gel(powers, 1) = identity;
    for (long exponent = 1; exponent < order; ++exponent)
        gel(powers, exponent + 1) =
            qfbcomp_i(gel(powers, exponent), generator);
    return powers;
}

static void
distance_encoding(GEN distance, double *mantissa, long *exponent)
{
    if (!signe(distance))
    {
        *mantissa = 0.0;
        *exponent = 0;
        return;
    }
    *exponent = gexpo(distance);
    *mantissa = gtodouble(gmul2n(distance, -*exponent));
}

static void
print_distance(FILE *output, GEN distance)
{
    double mantissa;
    long exponent;
    distance_encoding(distance, &mantissa, &exponent);
    fprintf(output, "[%.17g,%ld]", mantissa, exponent);
}

static GEN
round_polynomial(
    GEN approximate, long *global_exponent,
    GEN *maximum_real_distance, GEN *maximum_imaginary_distance,
    FILE *metrics, const char *label)
{
    long degree = degpol(approximate);
    GEN exact;
    *maximum_real_distance = gen_0;
    *maximum_imaginary_distance = gen_0;

    fprintf(metrics, "ROUNDING character=%s", label);
    for (long coefficient = 0; coefficient <= degree; ++coefficient)
    {
        GEN value = gel(approximate, coefficient + 2);
        GEN real_part = real_i(value);
        GEN imaginary_part = imag_i(value);
        long ignored;
        GEN nearest = grndtoi(real_part, &ignored);
        GEN real_distance = gabs(gsub(real_part, nearest), DEFAULTPREC);
        GEN imaginary_distance = gabs(imaginary_part, DEFAULTPREC);
        if (gcmp(real_distance, *maximum_real_distance) > 0)
            *maximum_real_distance = real_distance;
        if (gcmp(imaginary_distance, *maximum_imaginary_distance) > 0)
            *maximum_imaginary_distance = imaginary_distance;
        fprintf(metrics, " c%ld.real=", coefficient);
        print_distance(metrics, real_distance);
        fprintf(metrics, " c%ld.imag=", coefficient);
        print_distance(metrics, imaginary_distance);
    }
    fputc('\n', metrics);
    fflush(metrics);
    exact = grndtoi(approximate, global_exponent);
    return exact;
}

int
main(int argc, char **argv)
{
    if (argc < 3 || argc > 5)
    {
        fprintf(
            stderr,
            "usage: %s <base-polynomial> <output.gp> [safety-bits] "
            "[character-index]\n"
            "  character-index 1..6 selects one; 7 selects 2..6\n",
            argv[0]);
        return 2;
    }

    long safety_bits = argc == 4 ? atol(argv[3]) : 768;
    if (argc == 5) safety_bits = atol(argv[3]);
    if (safety_bits < 256)
    {
        fprintf(stderr, "safety-bits must be at least 256\n");
        return 2;
    }
    long selected_character = argc == 5 ? atol(argv[4]) : 0;
    if (selected_character < 0 || selected_character > 7)
    {
        fprintf(
            stderr,
            "character-index must be 1..7, or omitted for all\n");
        return 2;
    }

    pari_init_opts(
        1L << 29, 1048576,
        INIT_JMPm | INIT_SIGm | INIT_DFTm | INIT_noIMTm);
    paristack_setsize(1L << 29, 1L << 33);
    setrand(gen_1);
    double total_started = monotonic_seconds();

    GEN base_polynomial = gp_read_str(argv[1]);
    GEN K = Buchall(base_polynomial, nf_FORCE, DEFAULTPREC);
    GEN D = nf_get_disc(bnf_get_nf(K));
    GEN p = stoi(5);
    if (nf_get_degree(bnf_get_nf(K)) != 2
        || nf_get_r2(bnf_get_nf(K)) != 1)
        pari_err_BUG("CM constructor: base field is not imaginary quadratic");
    int certified = bnfcertify0(K, 0) == 1;
    if (!certified) pari_err_BUG("CM constructor: bnfcertify failed");

    GEN quadratic_group = quadclassunit0(D, 0, NULL, DEFAULTPREC);
    GEN cyc = bnf_get_cyc(K);
    GEN base_generators = bnf_get_gen(K);
    GEN class_number = bnf_get_no(K);
    if (!gequal(gel(quadratic_group, 1), class_number)
        || !gequal(gel(quadratic_group, 2), cyc))
        pari_err_BUG("CM constructor: bnf/quadclassunit mismatch");
    if (lg(cyc) != 4 || !equaliu(p, 5))
        pari_err_BUG("CM constructor: expected class-group rank three");
    long orders[3];
    long h = 1;
    for (long i = 0; i < 3; ++i)
    {
        orders[i] = itos(gel(cyc, i + 1));
        if (orders[i] % 5) pari_err_BUG("CM constructor: order not divisible by 5");
        h *= orders[i];
    }
    if (!equaliu(class_number, (ulong)h))
        pari_err_BUG("CM constructor: class number does not fit in long");

    struct gpq_data eta;
    initialize_double_eta(D, &eta);
    finish_double_eta_initialization(&eta);
    double gain = 12.0 * (eta.p + 1.0) * (eta.q + 1.0)
        / ((eta.p - 1.0) * (eta.q - 1.0));

    double pi_sqrt_D = M_PI * sqrt(fabs(gtodouble(D)));
    double kernel_size = h / 5.0;
    /* Explicit working height: leading j term divided by the actual
     * invariant gain, a trace allowance, and 32 nats for non-leading terms.
     * Five roots and the binomial bound 2^5 give the coefficient bound.
     */
    double root_height_nats =
        pi_sqrt_D / gain + log(kernel_size) + 32.0;
    double coefficient_height_nats = 5.0 * root_height_nats + 5.0 * M_LN2;
    long height_bits = (long)ceil(coefficient_height_nats / M_LN2);
    long precision_bits = height_bits + safety_bits;
    long precision = nbits2prec(precision_bits);
    eta.sqd = sqrtr_abs(itor(D, precision));

    GEN identity = principal_form(D);
    GEN generators = cgetg(4, t_VEC);
    GEN powers = cgetg(4, t_VEC);
    for (long i = 1; i <= 3; ++i)
    {
        gel(generators, i) = form_from_base_ideal(gel(base_generators, i), D);
        gel(powers, i) =
            form_powers(gel(generators, i), orders[i - 1], identity);
    }

    GEN sums[6];
    long counts[6][5] = {{0}};
    for (long character = 0; character < 6; ++character)
    {
        if ((selected_character >= 1 && selected_character <= 6
             && character + 1 != selected_character)
            || (selected_character == 7 && character == 0))
        {
            sums[character] = NULL;
            continue;
        }
        sums[character] = cgetg(6, t_VEC);
        for (long coset = 1; coset <= 5; ++coset)
            gel(sums[character], coset) = gen_0;
    }

    FILE *metrics = stdout;
    fprintf(
        metrics,
        "CONFIG D=%ld h=%ld cyc=[%ld,%ld,%ld] p=%ld q=%ld "
        "gain=%.17g pi_sqrt_D=%.17g kernel=%ld "
        "root_height_nats=%.17g coefficient_height_nats=%.17g "
        "height_bits=%ld safety_bits=%ld precision_bits=%ld "
        "precision_decimal_digits=%ld certified=%d selected_character=%ld\n",
        itos(D), h, orders[0], orders[1], orders[2], eta.p, eta.q,
        gain, pi_sqrt_D, h / 5, root_height_nats,
        coefficient_height_nats, height_bits, safety_bits, precision_bits,
        (long)floor(precision_bits * log10(2.0)), certified,
        selected_character);
    pari_printf("FORM_GENERATORS=%Ps\n", generators);
    fflush(metrics);

    double evaluation_started = monotonic_seconds();
    long evaluated = 0;
    for (long e1 = 0; e1 < orders[0]; ++e1)
        for (long e2 = 0; e2 < orders[1]; ++e2)
            for (long e3 = 0; e3 < orders[2]; ++e3)
            {
                GEN form = qfbcomp_i(
                    qfbcomp_i(gmael(powers, 1, e1 + 1),
                              gmael(powers, 2, e2 + 1)),
                    gmael(powers, 3, e3 + 1));
                GEN value = double_eta_value(form, &eta);
                ++evaluated;
                for (long character = 0; character < 6; ++character)
                {
                    if (!sums[character]) continue;
                    long coset = (
                        CHARACTER_COORDS[character][0] * (e1 % 5)
                        + CHARACTER_COORDS[character][1] * (e2 % 5)
                        + CHARACTER_COORDS[character][2] * (e3 % 5)) % 5;
                    gel(sums[character], coset + 1) =
                        gadd(gel(sums[character], coset + 1), value);
                    ++counts[character][coset];
                }
                if (evaluated % 1000 == 0)
                {
                    fprintf(
                        metrics, "PROGRESS evaluated=%ld/%ld seconds=%.6f\n",
                        evaluated, h, monotonic_seconds() - evaluation_started);
                    fflush(metrics);
                }
            }
    double evaluation_seconds = monotonic_seconds() - evaluation_started;
    fprintf(
        metrics, "EVALUATION_DONE classes=%ld seconds=%.9f\n",
        evaluated, evaluation_seconds);

    GEN raw_polynomials[6], reduced_polynomials[6];
    long rounding_exponents[6];
    GEN maximum_real_distances[6], maximum_imaginary_distances[6];
    double extraction_seconds[6];
    long polynomial_variable = fetch_user_var("y");

    for (long character = 0; character < 6; ++character)
    {
        if (!sums[character]) continue;
        double started = monotonic_seconds();
        for (long coset = 0; coset < 5; ++coset)
            if (counts[character][coset] != h / 5)
                pari_err_BUG("CM constructor: wrong kernel coset size");
        GEN approximate = roots_to_pol(sums[character], polynomial_variable);
        raw_polynomials[character] = round_polynomial(
            approximate, &rounding_exponents[character],
            &maximum_real_distances[character],
            &maximum_imaginary_distances[character], metrics,
            CHARACTER_LABELS[character]);
        if (rounding_exponents[character] > -128)
            pari_err_BUG("CM constructor: insufficient rounding separation");
        if (degpol(raw_polynomials[character]) != 5
            || !gequal1(leading_coeff(raw_polynomials[character])))
            pari_err_BUG("CM constructor: trace polynomial not monic quintic");
        /* Exact cosmetic reduction, strictly after integer recognition. */
        reduced_polynomials[character] =
            polredbest(raw_polynomials[character], 0);
        extraction_seconds[character] = monotonic_seconds() - started;
        pari_printf(
            "CHARACTER_DONE label=%s round_exponent=%ld raw=%Ps reduced=%Ps "
            "seconds=%.9f\n",
            CHARACTER_LABELS[character], rounding_exponents[character],
            raw_polynomials[character], reduced_polynomials[character],
            extraction_seconds[character]);
    }

    FILE *output = fopen(argv[2], "wb");
    if (!output) pari_err_FILE("CM constructor output", argv[2]);
    fprintf(output, "[1,%d,5,", PARI_VERSION_CODE);
    pari_fprintf(output, "%Ps,%Ps,", base_polynomial, D);
    pari_fprintf(
        output, "[%Ps,%Ps,%Ps,%d],",
        cyc, class_number, base_generators, certified);
    fprintf(
        output, "[\"double_eta_schertz\",%ld,%ld,%.17g],",
        eta.p, eta.q, gain);
    fprintf(
        output,
        "[%.17g,%ld,%.17g,%.17g,%ld,%ld,%ld,%ld],\n[",
        pi_sqrt_D, h / 5, root_height_nats, coefficient_height_nats,
        height_bits, safety_bits, precision_bits,
        (long)floor(precision_bits * log10(2.0)));
    int first_output_entry = 1;
    for (long character = 0; character < 6; ++character)
    {
        if (!sums[character]) continue;
        if (!first_output_entry) fprintf(output, ",\n");
        first_output_entry = 0;
        fprintf(output, "[\"%s\",[%ld,%ld,%ld]~ ,",
                CHARACTER_LABELS[character],
                CHARACTER_COORDS[character][0],
                CHARACTER_COORDS[character][1],
                CHARACTER_COORDS[character][2]);
        pari_fprintf(
            output, "%Ps,%Ps,%ld,",
            raw_polynomials[character], reduced_polynomials[character],
            rounding_exponents[character]);
        print_distance(output, maximum_real_distances[character]);
        fputc(',', output);
        print_distance(output, maximum_imaginary_distances[character]);
        fprintf(
            output, ",[%ld,%ld,%ld,%ld,%ld],%.9f]",
            counts[character][0], counts[character][1],
            counts[character][2], counts[character][3],
            counts[character][4], extraction_seconds[character]);
    }
    fprintf(
        output, "],%.9f,%.9f]\n",
        evaluation_seconds, monotonic_seconds() - total_started);
    if (fclose(output) != 0) pari_err_FILE("CM constructor output", argv[2]);

    fprintf(
        metrics, "CM_CONSTRUCTION_COMPLETE output=%s evaluation_seconds=%.9f "
        "total_seconds=%.9f\n",
        argv[2], evaluation_seconds, monotonic_seconds() - total_started);
    pari_close();
    return 0;
}
