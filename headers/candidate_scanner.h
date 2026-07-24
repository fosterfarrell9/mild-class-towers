// MIT License

#ifndef CANDIDATE_SCANNER_H
#define CANDIDATE_SCANNER_H

/**
 * Run the standalone imaginary-quadratic candidate-scanner CLI.
 *
 * This inexpensive filter enumerates negative fundamental field
 * discriminants, computes quadratic class groups, and retains fields with the
 * requested p-class rank. It never invokes secondary norms or mildness tests.
 */
int my_candidate_scanner_main(int argc, char **argv);

/** Print reproducible `--example-result` inputs from a candidate GP record. */
int my_candidate_inputs_main(int argc, char **argv);

/** Run focused fundamental-discriminant and class-rank scanner regressions. */
int my_candidate_scanner_self_test(void);

#endif
