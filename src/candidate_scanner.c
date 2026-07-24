// MIT License

/**
 * @file candidate_scanner.c
 * @brief Resumable p-class-rank filter for imaginary quadratic fields.
 *
 * The scan parameter is `n=|D_K|`. A value is examined exactly when `-n` is a
 * negative fundamental field discriminant, as decided by PARI's exact
 * `unegisfundamental`. For each such field, `quadclassunit` supplies the class
 * number and cyclic invariants; the p-rank is the number of invariants
 * divisible by p. This module is candidate discovery only and has no
 * dependency on the audited secondary-norm pipeline.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pari/pari.h>
#include "../headers/candidate_scanner.h"

#define CANDIDATE_FORMAT_VERSION 1
#define DEFAULT_CHECKPOINT_EVERY 1000000
#define DEFAULT_PROGRESS_SECONDS 10

typedef struct {
    long p;
    long rank;
    long min_abs;
    long max_abs;
    long checkpoint_every;
    long progress_seconds;
    long stop_after;
    const char *output;
    int resume;
} CandidateOptions;

typedef struct {
    long last_abs;
    long examined;
    GEN *items;
    long count;
    long capacity;
} CandidateState;

static double
candidate_elapsed(struct timespec start, struct timespec end)
{
    return (double)(end.tv_sec - start.tv_sec)
        + (double)(end.tv_nsec - start.tv_nsec) / 1e9;
}

static void
candidate_print_progress(
    const CandidateOptions *options, const CandidateState *state,
    struct timespec started, long considered)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    double elapsed = candidate_elapsed(started, now);
    double total =
        (double)(options->max_abs - options->min_abs + 1);
    double completed = state->last_abs < options->min_abs
        ? 0.0
        : (double)(state->last_abs - options->min_abs + 1);
    long seconds = (long)elapsed;
    pari_printf(
        "SCAN_PROGRESS |D|=%ld percent=%.2f fundamental=%ld "
        "candidates=%ld elapsed=%02ld:%02ld:%02ld rate=%.2f_abs_disc/s\n",
        state->last_abs, 100.0 * completed / total,
        state->examined, state->count,
        seconds / 3600, (seconds / 60) % 60, seconds % 60,
        elapsed > 0 ? considered / elapsed : 0.0);
    fflush(stdout);
}

static GEN
candidate_field(const char *name, GEN value)
{
    return mkvec2(strtoGENstr(name), value);
}

static long
candidate_p_rank(GEN cyc, GEN p)
{
    long rank = 0;
    for (long i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p)) ++rank;
    return rank;
}

static char *
candidate_polynomial(long discriminant)
{
    long residue = discriminant % 4;
    if (residue < 0) residue += 4;
    if (residue == 1)
        return pari_sprintf(
            "s^2-s+%ld", (1 - discriminant) / 4);
    return pari_sprintf("s^2+%ld", -discriminant / 4);
}

static void
candidate_state_append(CandidateState *state, GEN item)
{
    if (state->count == state->capacity)
    {
        long next = state->capacity ? 2 * state->capacity : 16;
        GEN *items = realloc(state->items, (size_t)next * sizeof(*items));
        if (!items) pari_err(e_MEM);
        state->items = items;
        state->capacity = next;
    }
    state->items[state->count++] = gclone(item);
}

static GEN
candidate_record(
    const CandidateOptions *options, const CandidateState *state,
    const char *status)
{
    GEN candidates = cgetg(state->count + 1, t_VEC);
    for (long i = 0; i < state->count; ++i)
        gel(candidates, i + 1) = state->items[i];

    GEN record = cgetg(11, t_VEC);
    gel(record, 1) = candidate_field("format_version", gen_1);
    gel(record, 2) = candidate_field("status", strtoGENstr(status));
    gel(record, 3) = candidate_field("p", stoi(options->p));
    gel(record, 4) = candidate_field("requested_p_rank", stoi(options->rank));
    gel(record, 5) = candidate_field("min_abs_discriminant", stoi(options->min_abs));
    gel(record, 6) = candidate_field("max_abs_discriminant", stoi(options->max_abs));
    gel(record, 7) = candidate_field(
        "last_abs_discriminant_considered", stoi(state->last_abs));
    gel(record, 8) = candidate_field(
        "fundamental_discriminants_examined", stoi(state->examined));
    gel(record, 9) = candidate_field("candidates_found", stoi(state->count));
    gel(record, 10) = candidate_field("candidates", candidates);
    return record;
}

static void
candidate_atomic_write(
    const CandidateOptions *options, const CandidateState *state,
    const char *status)
{
    pari_sp av = avma;
    size_t length = strlen(options->output) + 5;
    char *temporary = malloc(length);
    if (!temporary) pari_err(e_MEM);
    snprintf(temporary, length, "%s.tmp", options->output);
    FILE *file = fopen(temporary, "wb");
    if (!file) pari_err_FILE("candidate checkpoint", temporary);
    pari_fprintf(file, "%Ps\n", candidate_record(options, state, status));
    int flush_failed = fflush(file) != 0;
    int close_failed = fclose(file) != 0;
    if (flush_failed || close_failed)
        pari_err_FILE("candidate checkpoint", temporary);
    if (rename(temporary, options->output) != 0)
        pari_err_FILE("candidate checkpoint rename", options->output);
    free(temporary);
    avma = av;
}

static GEN
candidate_read_file(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (!file) pari_err_FILE("candidate checkpoint", path);
    if (fseek(file, 0, SEEK_END) != 0)
        pari_err_FILE("candidate checkpoint seek", path);
    long size = ftell(file);
    if (size < 0 || fseek(file, 0, SEEK_SET) != 0)
        pari_err_FILE("candidate checkpoint size", path);
    char *text = malloc((size_t)size + 1);
    if (!text) pari_err(e_MEM);
    if (fread(text, 1, (size_t)size, file) != (size_t)size)
        pari_err_FILE("candidate checkpoint read", path);
    fclose(file);
    text[size] = '\0';
    GEN result = gp_read_str(text);
    free(text);
    return result;
}

static void
candidate_load_checkpoint(
    const CandidateOptions *options, CandidateState *state)
{
    GEN record = candidate_read_file(options->output);
    if (typ(record) != t_VEC || glength(record) != 10
        || !equaliu(gmael(record, 1, 2), CANDIDATE_FORMAT_VERSION))
        pari_err(e_MISC, "candidate checkpoint schema mismatch");
    if (!equaliu(gmael(record, 3, 2), options->p)
        || !equaliu(gmael(record, 4, 2), options->rank)
        || !equaliu(gmael(record, 5, 2), options->min_abs)
        || !equaliu(gmael(record, 6, 2), options->max_abs))
        pari_err(e_MISC, "candidate checkpoint scan parameters differ");

    state->last_abs = itos(gmael(record, 7, 2));
    state->examined = itos(gmael(record, 8, 2));
    GEN candidates = gmael(record, 10, 2);
    if (typ(candidates) != t_VEC
        || glength(candidates) != itos(gmael(record, 9, 2)))
        pari_err(e_MISC, "candidate checkpoint count mismatch");
    for (long i = 1; i < lg(candidates); ++i)
        candidate_state_append(state, gel(candidates, i));
}

static long
candidate_parse_long(const char *option, const char *value)
{
    char *end = NULL;
    errno = 0;
    long parsed = strtol(value, &end, 10);
    if (errno || !value[0] || !end || *end || parsed <= 0)
    {
        fprintf(stderr, "%s requires a positive integer\n", option);
        return -1;
    }
    return parsed;
}

static int
candidate_parse_options(
    int argc, char **argv, CandidateOptions *options)
{
    memset(options, 0, sizeof(*options));
    options->checkpoint_every = DEFAULT_CHECKPOINT_EVERY;
    options->progress_seconds = DEFAULT_PROGRESS_SECONDS;
    for (int i = 2; i < argc; ++i)
    {
        if (strcmp(argv[i], "--resume") == 0)
            options->resume = 1;
        else if (i + 1 < argc)
        {
            long value;
            if (strcmp(argv[i], "--output") == 0)
                options->output = argv[++i];
            else if (strcmp(argv[i], "--prime") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->p = value;
            }
            else if (strcmp(argv[i], "--rank") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->rank = value;
            }
            else if (strcmp(argv[i], "--min-abs-disc") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->min_abs = value;
            }
            else if (strcmp(argv[i], "--max-abs-disc") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->max_abs = value;
            }
            else if (strcmp(argv[i], "--checkpoint-every") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->checkpoint_every = value;
            }
            else if (strcmp(argv[i], "--progress-seconds") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->progress_seconds = value;
            }
            else if (strcmp(argv[i], "--stop-after") == 0)
            {
                value = candidate_parse_long(argv[i], argv[i + 1]);
                ++i;
                if (value < 0) return 0;
                options->stop_after = value;
            }
            else
            {
                fprintf(stderr, "Unknown candidate option: %s\n", argv[i]);
                return 0;
            }
        }
        else
        {
            fprintf(stderr, "Missing candidate option value: %s\n", argv[i]);
            return 0;
        }
    }
    if (!options->output || options->p < 3 || !(options->p & 1)
        || options->rank <= 0 || options->min_abs <= 0
        || options->max_abs < options->min_abs)
    {
        fprintf(stderr,
                "Usage: %s --scan-candidates --prime <odd-p> --rank <r> "
                "--min-abs-disc <n> --max-abs-disc <n> --output <file.gp> "
                "[--checkpoint-every <n>] [--progress-seconds <n>] "
                "[--resume] [--stop-after <n>]\n",
                argv[0]);
        return 0;
    }
    return 1;
}

int
my_candidate_scanner_main(int argc, char **argv)
{
    CandidateOptions options;
    if (!candidate_parse_options(argc, argv, &options))
        return EXIT_FAILURE;
    pari_init(1L << 28, 1048576);
    GEN p = stoi(options.p);
    if (!uisprime((ulong)options.p))
    {
        fprintf(stderr, "--prime must be prime\n");
        pari_close();
        return EXIT_FAILURE;
    }

    CandidateState state = {0};
    state.last_abs = options.min_abs - 1;
    if (options.resume)
    {
        candidate_load_checkpoint(&options, &state);
        pari_printf(
            "RESUME last=%ld fundamental=%ld candidates=%ld next=%ld\n",
            state.last_abs, state.examined, state.count, state.last_abs + 1);
    }
    long start_abs = state.last_abs + 1;
    long considered = 0;
    struct timespec started;
    clock_gettime(CLOCK_MONOTONIC, &started);
    struct timespec next_progress = started;
    next_progress.tv_sec += options.progress_seconds;
    pari_printf(
        "SCAN_START p=%ld rank=%ld min_abs=%ld max_abs=%ld "
        "checkpoint_every=%ld progress_seconds=%ld PARI=%Ps\n",
        options.p, options.rank, options.min_abs, options.max_abs,
        options.checkpoint_every, options.progress_seconds, pari_version());
    candidate_print_progress(&options, &state, started, considered);

    for (long n = start_abs; n <= options.max_abs; ++n)
    {
        pari_sp av = avma;
        state.last_abs = n;
        ++considered;
        if (unegisfundamental((ulong)n))
        {
            ++state.examined;
            GEN data = quadclassunit0(stoi(-n), 0, NULL, DEFAULTPREC);
            GEN cyc = gel(data, 2);
            long rank = candidate_p_rank(cyc, p);
            if (rank == options.rank)
            {
                char *polynomial = candidate_polynomial(-n);
                GEN item = cgetg(7, t_VEC);
                gel(item, 1) = stoi(-n);
                gel(item, 2) = p;
                gel(item, 3) = cyc;
                gel(item, 4) = gel(data, 1);
                gel(item, 5) = stoi(rank);
                gel(item, 6) = strtoGENstr(polynomial);
                candidate_state_append(&state, item);
                pari_free(polynomial);
                pari_printf(
                    "CANDIDATE D=%ld class_group=%Ps h=%Ps p_rank=%ld\n",
                    -n, cyc, gel(data, 1), rank);
                fflush(stdout);
            }
        }
        avma = av;

        if (considered % options.checkpoint_every == 0)
        {
            candidate_atomic_write(&options, &state, "IN_PROGRESS");
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = candidate_elapsed(started, now);
            pari_printf(
                "CHECKPOINT |D|=%ld fundamental=%ld candidates=%ld "
                "elapsed=%.3f\n",
                state.last_abs, state.examined, state.count, elapsed);
            fflush(stdout);
        }
        if (considered % 10000 == 0)
        {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            if (now.tv_sec > next_progress.tv_sec
                || (now.tv_sec == next_progress.tv_sec
                    && now.tv_nsec >= next_progress.tv_nsec))
            {
                candidate_print_progress(
                    &options, &state, started, considered);
                next_progress = now;
                next_progress.tv_sec += options.progress_seconds;
            }
        }
        if (options.stop_after && considered >= options.stop_after
            && n < options.max_abs)
        {
            candidate_atomic_write(&options, &state, "IN_PROGRESS");
            pari_printf("SCAN PAUSED last=%ld\n", state.last_abs);
            goto finished;
        }
    }
    candidate_atomic_write(&options, &state, "COMPLETE");
    candidate_print_progress(&options, &state, started, considered);

finished:
    {
        struct timespec ended;
        clock_gettime(CLOCK_MONOTONIC, &ended);
        double elapsed = candidate_elapsed(started, ended);
        pari_printf(
            "SCAN SUMMARY last=%ld fundamental=%ld candidates=%ld "
            "elapsed=%.3f rate=%.2f_abs_disc/s\n",
            state.last_abs, state.examined, state.count, elapsed,
            elapsed > 0 ? considered / elapsed : 0.0);
    }
    for (long i = 0; i < state.count; ++i) gunclone(state.items[i]);
    free(state.items);
    pari_close();
    return EXIT_SUCCESS;
}

int
my_candidate_inputs_main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s --candidate-inputs <candidates.gp>\n", argv[0]);
        return EXIT_FAILURE;
    }
    pari_init(1L << 26, 1048576);
    GEN record = candidate_read_file(argv[2]);
    if (typ(record) != t_VEC || glength(record) != 10)
        pari_err(e_MISC, "candidate input file schema mismatch");
    GEN candidates = gmael(record, 10, 2);
    for (long i = 1; i < lg(candidates); ++i)
        pari_printf(
            "p=%Ps D=%Ps polynomial='%s' command: "
            "./build/massey --example-result RESULT.gp %Ps '%s'\n",
            gmael(candidates, i, 2), gmael(candidates, i, 1),
            GSTR(gmael(candidates, i, 6)),
            gmael(candidates, i, 2), GSTR(gmael(candidates, i, 6)));
    pari_close();
    return EXIT_SUCCESS;
}

int
my_candidate_scanner_self_test(void)
{
    pari_init(1L << 27, 1048576);
    const long fundamental[] = {3, 4, 7, 8, 11, 11203620, 18397407};
    const long nonfundamental[] = {1, 5, 9, 12, 16, 11203621};
    for (size_t i = 0; i < sizeof(fundamental) / sizeof(*fundamental); ++i)
        if (!unegisfundamental((ulong)fundamental[i]))
            pari_err(e_MISC, "fundamental-discriminant recognition failed");
    for (size_t i = 0; i < sizeof(nonfundamental) / sizeof(*nonfundamental); ++i)
        if (unegisfundamental((ulong)nonfundamental[i]))
            pari_err(e_MISC, "nonfundamental discriminant accepted");

    GEN p = stoi(5);
    GEN first = quadclassunit0(stoi(-11203620), 0, NULL, DEFAULTPREC);
    GEN second = quadclassunit0(stoi(-18397407), 0, NULL, DEFAULTPREC);
    GEN rejected = quadclassunit0(stoi(-3), 0, NULL, DEFAULTPREC);
    if (!gequal(gel(first, 2), mkvec3(stoi(10), stoi(10), stoi(10)))
        || candidate_p_rank(gel(first, 2), p) != 3
        || candidate_p_rank(gel(second, 2), p) != 3
        || candidate_p_rank(gel(rejected, 2), p) == 3)
        pari_err(e_MISC, "candidate class-rank regression failed");
    pari_printf(
        "CANDIDATE_SELF_TEST PASS "
        "D=-11203620 Cl=%Ps; D=-18397407 Cl=%Ps; D=-3 rejected\n",
        gel(first, 2), gel(second, 2));
    pari_close();
    return EXIT_SUCCESS;
}
