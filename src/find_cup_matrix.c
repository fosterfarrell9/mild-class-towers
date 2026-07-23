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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pari/pari.h>
#include "../headers/find_cup_matrix.h"
#include "../headers/artin_symbol.h"
#include "../headers/misc_functions.h"
#include "../headers/secondary_norm.h"

// Debug macros
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define MY_DEBUGLEVEL 0
#define DEBUG_PRINT(level, ...) \
    do { if (MY_DEBUGLEVEL >= (level)) pari_printf(__VA_ARGS__); } while (0)

static int
my_diagnostics_enabled(void)
{
    const char *value = getenv("MASSEY_DIAGNOSTICS");
    return value && strcmp(value, "1") == 0;
}

static void
my_diagnostic_fail(const char *message)
{
    fprintf(stderr, "MASSEY_DIAGNOSTICS failure: %s\n", message);
    exit(EXIT_FAILURE);
}

/* [c, A, lambda, u, x_H90, c_over_x, D_scale] */
static GEN
my_extension_diagnostic(GEN K_ext, GEN K, GEN p, GEN H, long k)
{
    pari_sp av = avma;
    GEN Labs = gmael(K_ext, k, 1);
    GEN Lrel = gmael(K_ext, k, 2);
    GEN sigma_H90 = gmael(K_ext, k, 3);

    /* Complete Lrel with its absolute nf before calling rnfcycaut. */
    (void)rnfidealup0(Lrel, idealhnf0(K, gen_1, NULL), 1);

    GEN c = my_subgroup_character(K, H, p);
    GEN A = zerocol(glength(c));
    GEN cyc = bnf_get_cyc(K), generators = bnf_get_gen(K);
    long i, j = 1, first = 0, p_int = itos(p);

    for (i = 1; i < lg(cyc); ++i)
        if (dvdii(gel(cyc, i), p))
        {
            long value = my_Artin_symbol(
                Labs, Lrel, K, gel(generators, i), p_int);
            gel(A, j) = stoi((value % p_int + p_int) % p_int);
            if (!first && signe(gel(c, j))) first = j;
            ++j;
        }

    if (!first)
        my_diagnostic_fail("normalized subgroup character is zero");
    GEN lambda = Fp_div(gel(A, first), gel(c, first), p);
    if (!signe(lambda))
        my_diagnostic_fail("Artin character is zero");
    for (i = 1; i < lg(c); ++i)
        if (!equalii(gel(A, i), Fp_mul(lambda, gel(c, i), p)))
            my_diagnostic_fail("Artin character is not a scalar multiple of subgroup character");

    long u_long = my_sigma_exponent(Labs, Lrel, sigma_H90, p);
    GEN u = stoi(u_long);
    GEN u_inverse = Fp_inv(u, p);
    GEN x_H90 = cgetg(lg(A), t_COL);
    for (i = 1; i < lg(A); ++i)
        gel(x_H90, i) = Fp_mul(u_inverse, gel(A, i), p);

    GEN c_over_x = Fp_div(u, lambda, p);
    GEN D_scale = Fp_sqr(c_over_x, p);
    GEN result = mkvecn(7, c, A, lambda, u, x_H90, c_over_x, D_scale);

    pari_printf("\nMASSEY_DIAGNOSTICS extension k=%ld\n", k);
    pari_printf("  subgroup H_k = %Ps\n", H);
    pari_printf("  subgroup character line c_k = %Ps\n", gtovec(c));
    pari_printf("  Artin character A_k = %Ps\n", gtovec(A));
    pari_printf("  lambda_k = %Ps\n", lambda);
    pari_printf("  sigma exponent u_k = %Ps\n", u);
    pari_printf("  H90 character x_H90,k = %Ps\n", gtovec(x_H90));
    pari_printf("  scale c_k / x_H90,k = %Ps\n", c_over_x);
    pari_printf("  D scaling factor = %Ps\n", D_scale);
    return gerepilecopy(av, result);
}

//-----------------------------------------------------------------------------

//Defines a matrix over F_2 with index (i*k, j) corresponding to 
//< x_i\cup x_k, (a_j, J_j)>

//------------------------------------------------
// Parallel version of my_relations
//------------------------------------------------

/*
GP;install("compute_my_relations","Gp","./stdin.so");
*/

// Wrapper function for parallel computation
GEN compute_my_relations(long i, GEN args) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: compute_my_relations\n--------------------------\n\n");
    pari_sp av = avma;

    GEN K_ext   = gel(args, 1);
    GEN K       = gel(args, 2);
    GEN p       = gel(args, 3); 
    GEN gp_rk   = gel(args, 4);
    GEN Ja_vect = gel(args, 5);
    GEN gr_rk   = gel(args, 6);

    
    GEN result = zerovec(2); // A vector to store results
    
    int p_rk = itos(gp_rk), r_rk = itos(gr_rk);
    
    
    GEN s_cup_matrix = zerovec(r_rk);
    GEN s_cup_matrix_full = zerovec(r_rk);

    GEN Labs_cup, Lrel_cup, sigma_cup, I_vect, I_rel, NIpJ, Labs, Lrel, Lbnr_cup;
    long j, k;

    for (j = 1; j <= r_rk; ++j) {
        gel(s_cup_matrix, j) = zerovec((p_rk * (p_rk + 1)) / 2);
        gel(s_cup_matrix_full, j) = zerovec(p_rk * p_rk);
    }

    Labs_cup = gmael(K_ext, i, 1);
    Lrel_cup = gmael(K_ext, i, 2);
    sigma_cup = gmael(K_ext, i, 3);
    Lbnr_cup = gmael(K_ext, i, 4);
    
    //--------------------------------------------------------------------------------
    // Artin symbol test 
    //my_test_artin_symbol(Labs_cup, Lrel_cup, K, itos(p), sigma_cup);
    //--------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------
    // Find the I's
    I_vect = my_H90_vect_2(Labs_cup, Lrel_cup, Lbnr_cup, K, sigma_cup, Ja_vect, p, 1);
    //--------------------------------------------------------------------------------

    for (j = 1; j < r_rk + 1; ++j) {
        
        I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
        if (itos(p) == 2) {
            NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
        } else {
            NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
        }
        for (k = 1; k < p_rk + 1; ++k) {
            
            Labs = gel(gel(K_ext, k), 1);
            Lrel = gel(gel(K_ext, k), 2);
            
            if (i < k) {
                
                gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)) = stoi(my_Artin_symbol(Labs, Lrel, K, idealred0(K, NIpJ, NULL), itos(p))%itos(p));
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1));
                
            }
            if (i == k) {
                // The Bockstein case ...
                gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)) = stoi(my_Artin_symbol(Labs, Lrel, K, gel(gel(Ja_vect, j), 2), itos(p))%itos(p));
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = gmael2(s_cup_matrix, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1));
                
            }
            if (i > k) {
                
                gmael2(s_cup_matrix_full, j, p_rk * (k - 1) + i) = stoi((itos(p)-my_Artin_symbol(Labs, Lrel, K, idealred0(K, NIpJ, NULL), itos(p)))%itos(p));
                
            }
        }
    }

    gel(result, 1) = s_cup_matrix;
    
    gel(result, 2) = s_cup_matrix_full;
    
    DEBUG_PRINT(1, "\n--------------------------\nEnd: compute_my_relations[i=%d<=%d]\n--------------------------\n\n", i, p_rk);
    result = gerepilecopy(av, result);
    return result;
}

int my_relations_par(GEN K_ext, GEN K, GEN p, int p_rk, GEN Ja_vect, int r_rk) {
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_relations_par\n--------------------------\n\n");
    //GEN NIpJ, I_vect, I_rel, Labs, Lrel, sigma, Labs_cup, Lrel_cup, sigma_cup;

    int nr_col = (p_rk * (p_rk + 1) / 2), nr_col_full = p_rk * p_rk;
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row), cup_matrix_full = zerovec(nr_row);
    
    for (int j = 1; j < nr_row + 1; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
        gel(cup_matrix_full, j) = zerovec(nr_col_full);
    }

    long i, taskid, pending;
    GEN done, res, res_full;

    struct pari_mt pt;
    
    GEN args = mkvecn(6, K_ext, K, p, stoi(p_rk), Ja_vect, stoi(r_rk));
    //DEBUG_PRINT(1, "strtoclosure: %Ps\n", strtoclosure("_worker", 1, args));
    mt_queue_start(&pt, strtoclosure("_worker", 1, args));
    
    for (i = 1; i <= p_rk || pending; i++)
    { 
        DEBUG_PRINT(1, "for i = %ld <= %d\n", i, p_rk);
        mt_queue_submit(&pt, i, i<=p_rk? mkvecs(i): NULL);
        done = mt_queue_get(&pt, &taskid, &pending);
        if (done) {
            DEBUG_PRINT(1, "Done i = %ld \n", i);
            res = gel(done, 1);
            res_full = gel(done, 2);

            for (int j = 1; j <= r_rk; ++j) {
                for (int k = 1; k <= nr_col; ++k) {
                    gmael2(cup_matrix, j, k) = gadd(gmael2(cup_matrix, j, k), gmael2(res, j, k));
                }
                for (int k = 1; k <= nr_col_full; ++k) {
                    gmael2(cup_matrix_full, j, k) = gadd(gmael2(cup_matrix_full, j, k), gmael2(res_full, j, k));
                }
            }
        }
        
    }
    
    mt_queue_end(&pt); /* end parallelism */


    DEBUG_PRINT(1, ANSI_COLOR_YELLOW "\n\nFull Cup Matrix:  (for (i,k) with i>k we have - (x_i cup x_k) instead of x_i cup x_k in order to see the symmetry).\n\n" ANSI_COLOR_RESET);
    for (int j = 1; j < nr_row + 1; ++j) {
        DEBUG_PRINT(1, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix_full, j));
    }

    int sym_check = 1;
    GEN fail;
    for (i = 1; i <= p_rk; i++) {
        for (int j = 1; j < r_rk; j++) {
            for (int k = 1; k <= p_rk; k++) {
                if (!gequal(gmael2(cup_matrix_full, j, p_rk * (k - 1) + i), gmael2(cup_matrix_full, j, p_rk * (i - 1) + k))) {
                    fail = mkvec2s(i, k);
                    sym_check = 0;
                    break;
                }
            }
        }
    }
    switch (sym_check) {
    case 1:
        DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nSymmetry test passed\n\n" ANSI_COLOR_RESET);
        break;
    case 0:
        DEBUG_PRINT(1, ANSI_COLOR_RED "\nSymmetry test failed at: %Ps\n\n" ANSI_COLOR_RESET, fail);
        pari_close();
        exit(111);
        break;
    default:
        break;
    }
    
    DEBUG_PRINT(0, "---------------------------------------------------------------------------------------------------------\n\n");
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (int j = 1; j <= nr_row; ++j) {
        DEBUG_PRINT(0, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    DEBUG_PRINT(0, ANSI_COLOR_RED "For indices (i,i) we have B(x_i) instead of x_i cup x_i to get the correct presentation of Q_2.\n\n" ANSI_COLOR_RESET);
    DEBUG_PRINT(0, ANSI_COLOR_CYAN "This determines the second quotient Q_2 for the lower p-central series and not only the Zassenhaus quotient ZQ_2.\n\n" ANSI_COLOR_RESET);
    DEBUG_PRINT(0, "---------------------------------------------------------------------------------------------------------\n\n");
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    DEBUG_PRINT(0, ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))), p);
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (int i = 1; i < glength(cup_hnf) + 1; ++i) {
            DEBUG_PRINT(0, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        }
        DEBUG_PRINT(0, "\n");
        
        char letters[] = "abcdefghijklmnopqr";
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        DEBUG_PRINT(0, "\"\"\"[");
        for (int j = 1; j <= hnf_r_rk; j++) {
            for (int i = 1; i <= p_rk; ++i) {
                for (int k = i; k <= p_rk; k++) {
                    if (!gequal0(gel(gel(cup_hnf, j), (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1)))) {
                        if (i == k) {
                            DEBUG_PRINT(0, "%c^%ld", letters[i - 1], itos(p) * smodis(gneg(gmael2(cup_hnf, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1))), itos(p)));
                        }
                        if (i < k) {
                            DEBUG_PRINT(0, "%c_%c^%ld", letters[i - 1], letters[k - 1], smodis(gneg(gmael2(cup_hnf, j, (2 * p_rk - (i - 2)) * (i - 1) / 2 + k - (i - 1))), itos(p)));
                        }
                    }
                }
            }
            if (j == hnf_r_rk) {
                DEBUG_PRINT(0, "]\"\"\"\n");
            } else {
                DEBUG_PRINT(0, ", ");
            }
        }
    }
    else {
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        DEBUG_PRINT(0, "\"\"\"[]\"\"\"\n");
    }
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_relations_par\n--------------------------\n\n");
    DEBUG_PRINT(0, "---------------------------------------------------------------------------------------------------------\n\n");

    return mat_rk;
}

int my_relations (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk)
{
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_relations\n--------------------------\n\n");
    GEN NIpJ, I_vect, I_rel, Labs, Lrel, Labs_cup, Lrel_cup, Lbnr_cup, sigma_cup;

    //-------------------------------------------------
    // Define a matrix with r_rk rows and columns indexed by (i, k) for 1 <= i <= k <= p_rk. 

    int nr_col = (p_rk*(p_rk+1)/2), nr_col_full = p_rk*p_rk;
    int nr_row = r_rk;
    GEN cup_matrix = zerovec(nr_row), cup_matrix_full = zerovec(nr_row);
    
    int i, j, k;
    for (j=1; j<=nr_row; ++j) {
        gel(cup_matrix, j) = zerovec(nr_col);
        gel(cup_matrix_full, j) = zerovec(nr_col_full);
    }
    //-------------------------------------------------

    // <x_i cup x_k, (a_j, J_j)> = < x_k , ... > ( = Artin symbol) ----->  i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<=p_rk; ++i) {
        //DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "-----------\n\n\nStarting i = %d\n\n\n-------------\n" ANSI_COLOR_RESET, i);
        Labs_cup = gmael(K_ext, i, 1);
        Lrel_cup = gmael(K_ext, i, 2);
        sigma_cup = gmael(K_ext, i, 3);
        Lbnr_cup = gmael(K_ext, i, 4);

        // Artin symbol test
        //my_test_artin_symbol(Labs_cup, Lrel_cup, K, p_int, sigma_cup);
        // my_test_artin_on_norms(Labs_cup, Lrel_cup, K, p_int, sigma_cup);
        // I_vect corresp. to i:th extension
        // I_vect = my_find_I_vect_full(Labs_cup, Lrel_cup, K, sigma_cup, Ja_vect, p_int);

        I_vect = my_H90_vect_2(Labs_cup, Lrel_cup, Lbnr_cup, K, sigma_cup, Ja_vect, stoi(p_int), 1);
        
        for (j=1; j<=r_rk; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            // DEBUG_PRINT(1, ANSI_COLOR_MAGENTA "-----------\n\n\nStarting [i,j] = [%d, %d]\n\n\n-------------\n" ANSI_COLOR_RESET,i, j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_vect, j));
            // DEBUG_PRINT(1, "I, %d to rel\n\n", j);
            
            if (p_int == 2) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j), 2), rnfidealnormrel(Lrel_cup, I_rel));
            }
            else {
                NIpJ = rnfidealnormrel(Lrel_cup, I_rel);
            }
            
            for (k=1; k<=p_rk; ++k) {
                // DEBUG_PRINT(1, ANSI_COLOR_YELLOW "----------------\n\n\n\nStarting [i,j,k] = [%d, %d, %d]\n\n\n\n----------------\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                
                if (i<k) {
                    gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int)%p_int);
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1));
                }
                if (i==k) {
                    gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1)) = stoi(my_Artin_symbol(Labs, Lrel, K, gel(gel(Ja_vect, j), 2), p_int)%p_int);
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = gmael2(cup_matrix, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1));
                }
                if (i>k) {
                    gmael2(cup_matrix_full, j, p_rk*(k-1)+i) = stoi((p_int-my_Artin_symbol(Labs, Lrel, K, idealred0(K,NIpJ, NULL), p_int))%p_int);
                }
                
            }
            
        }
        
    }
    // ------------------------
    //Check symmetry:
    DEBUG_PRINT(1, ANSI_COLOR_YELLOW "Full Cup Matrix:  (for (i,k) with i>k we have - (x_i cup x_k) instead of x_i cup x_k in order to see the symmetry).\n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
        DEBUG_PRINT(1, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix_full, j));
    }
    int sym_check = 1;
    GEN fail;
    for (i = 1; i <= p_rk; i++)
    {
        for (j = 1; j < r_rk; j++)
        {
            for (k = 1; k <= p_rk; k++)
            {
                if (!gequal(gmael2(cup_matrix_full, j, p_rk*(k-1)+i), gmael2(cup_matrix_full, j, p_rk*(i-1)+k)))
                {
                    // gel(fail, 1) = stoi(i);
                    // gel(fail, 2) = stoi(k);
                    fail = mkvec2s(i, k);
                    
                    sym_check = 0;
                    break;
                }
                
            }
        }
    }
    switch (sym_check)
    {
    case 1:
        DEBUG_PRINT(1, ANSI_COLOR_GREEN "\nSymmetry test passed\n\n" ANSI_COLOR_RESET);
        break;
    case 0:
        DEBUG_PRINT(1, ANSI_COLOR_RED "\nSymmetry test failed at: %Ps\n\n" ANSI_COLOR_RESET, fail);
        pari_close();
        exit(111);
        break;
    default:
        break;
    }
    
    DEBUG_PRINT(0, "---------------------------------------------------------------------------------------------------------\n\n");
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<=nr_row; ++j) {
        DEBUG_PRINT(1, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_matrix, j));
    }
    DEBUG_PRINT(1, ANSI_COLOR_RED "For indices (i,i) we have B(x_i) instead of x_i cup x_i to get the correct presentation of Q_2 (Koch 7.23 vs Koch 7.24).\n\n" ANSI_COLOR_RESET);
    
    DEBUG_PRINT(0, ANSI_COLOR_CYAN "This determines the second quotient Q_2 for the lower p-central series and not only the Zassenhaus quotient ZQ_2.\n\n" ANSI_COLOR_RESET);
    DEBUG_PRINT(0, "---------------------------------------------------------------------------------------------------------\n\n");
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(cup_matrix)), p);
    DEBUG_PRINT(0, ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);

    if (mat_rk > 0) {
        GEN cup_hnf = FpM_red(hnf((ZM_copy(cup_matrix))),p);
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (i=1;i<glength(cup_hnf)+1;++i) {
            DEBUG_PRINT(0, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(cup_hnf, i));
        }
        
        DEBUG_PRINT(0, "\n\n");
        
        char letters[] = "abcdefghijklmnopqr";
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        int hnf_r_rk = glength(cup_hnf);
        DEBUG_PRINT(0, "\"\"\"[");
        for (j=1; j<=hnf_r_rk; j++) {
            for (i=1; i<=p_rk; ++i) {
                for (k=i; k<=p_rk; k++) {
                    if (!gequal0(gel(gel(cup_hnf, j), (2*p_rk-(i-2))*(i-1)/2+k-(i-1)))) {
                        if (i==k) { // note here that we use gneg to get the correct sign as in Koch 7.24
                            DEBUG_PRINT(0, "%c^%ld", letters[i-1], p_int*smodis(gneg(gmael2(cup_hnf, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1))), p_int));
                        }
                        if (i<k) { // note here that we use gneg to get the correct sign as in Koch 7.24
                            DEBUG_PRINT(0, "%c_%c^%ld", letters[i-1],letters[k-1], smodis(gneg(gmael2(cup_hnf, j, (2*p_rk-(i-2))*(i-1)/2+k-(i-1))), p_int));
                        }
                    }
                }
            }
            if (j==hnf_r_rk) {
                DEBUG_PRINT(0, "]\"\"\"\n");
            }
            else {
                DEBUG_PRINT(0, ", ");
            }
        }
    }
    else {
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Cup relations:  \n\n" ANSI_COLOR_RESET);
        DEBUG_PRINT(0, "\"\"\"[]\"\"\"\n");
    }
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_relations\n--------------------------\n\n");
    DEBUG_PRINT(1, "---------------------------------------------------------------------------------------------------------\n\n");
    
    return mat_rk;
}



int my_massey_matrix (GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups, int n)
{
    GEN NI, NIpJ, I_rel, Labs, Lrel, Labs_cup, Lrel_cup, Lbnr_cup, sigma_cup, I_prime_vect;
    int nr_col = p_rk*p_rk;
    int nr_row = r_rk;
    GEN massey_matrix = zerovec(nr_row);
    int diagnostics = my_diagnostics_enabled() && p_int > 3 && n == 2;
    int scaling_test = diagnostics && p_int == 5;
    GEN diagnostic_data = NULL;
    
    int i, j, k;
    for (j=1; j<(nr_row+1); ++j) {
        gel(massey_matrix, j) = zerovec(nr_col);
    }

    if (diagnostics)
    {
        if (glength(best_subgroups) != p_rk)
            my_diagnostic_fail("selected subgroup count does not equal p-rank");
        diagnostic_data = zerovec(p_rk);
        pari_printf("\nMASSEY_DIAGNOSTICS triple-product normalization\n");
        for (i = 1; i <= p_rk; ++i)
            gel(diagnostic_data, i) =
                my_extension_diagnostic(K_ext, K, p,
                                        gel(best_subgroups, i), i);
    }

    //DEBUG_PRINT(1, "cup_mat: %Ps\n\n", cup_matrix);
    // i:th extension, j:th ideal J, evaluate on k:th extension 
    for (i=1; i<p_rk+1; ++i) {
        DEBUG_PRINT(1, "Start round %d/%d\n\n", i, p_rk);
        Labs_cup = gel(gel(K_ext, i), 1);
        Lrel_cup = gel(gel(K_ext, i), 2);
        sigma_cup = gel(gel(K_ext, i), 3);
        Lbnr_cup = gel(gel(K_ext, i), 4);

        //--------------------------------------------------------------------------------
        // Find the I's
        I_prime_vect = my_H90_vect_2(Labs_cup, Lrel_cup, Lbnr_cup, K, sigma_cup, Ja_vect, p, n);
        DEBUG_PRINT(1, "I'_vect found\n\n");
        //--------------------------------------------------------------------------------

        GEN sigma_squared = NULL, I_prime_vect_sigma2 = NULL;
        GEN x_sigma2 = NULL, expected_factor = NULL;
        if (scaling_test)
        {
            sigma_squared = my_automorphism_power_checked(
                Labs_cup, Lrel_cup, K, sigma_cup, 2, p);

            GEN A_i = gmael(diagnostic_data, i, 2);
            GEN x_sigma = gmael(diagnostic_data, i, 5);
            GEN u_i = gmael(diagnostic_data, i, 4);
            GEN two = stoi(2);
            GEN inverse_two = Fp_inv(two, p);
            long u_sigma2_long =
                my_sigma_exponent(Labs_cup, Lrel_cup,
                                  sigma_squared, p);
            GEN u_sigma2 = stoi(u_sigma2_long);
            GEN expected_u_sigma2 = Fp_mul(two, u_i, p);
            if (!equalii(u_sigma2, expected_u_sigma2))
                my_diagnostic_fail(
                    "sigma squared has the wrong exponent relative to sigma_Artin");

            x_sigma2 = cgetg(lg(A_i), t_COL);
            for (long coordinate = 1;
                 coordinate < lg(A_i); ++coordinate)
                gel(x_sigma2, coordinate) =
                    Fp_mul(Fp_inv(u_sigma2, p),
                           gel(A_i, coordinate), p);

            for (long coordinate = 1;
                 coordinate < lg(x_sigma); ++coordinate)
                if (!equalii(
                        gel(x_sigma2, coordinate),
                        Fp_mul(inverse_two,
                               gel(x_sigma, coordinate), p)))
                    my_diagnostic_fail(
                        "sigma squared character does not scale by 2^(-1)");

            expected_factor = Fp_sqr(inverse_two, p);
            if (!equaliu(expected_factor, 4))
                my_diagnostic_fail(
                    "p=5 sigma-squared expected factor is not 4");

            pari_printf(
                "\nMASSEY_DIAGNOSTICS H90 scaling setup i=%d\n", i);
            pari_printf("  sigma power a = 2\n");
            pari_printf("  sigma^2 exponent relative to sigma_Artin = %Ps\n",
                        u_sigma2);
            pari_printf("  x_sigma = %Ps\n", gtovec(x_sigma));
            pari_printf("  x_sigma2 = %Ps\n", gtovec(x_sigma2));
            pari_printf("  character scale = %Ps\n", inverse_two);
            pari_printf("  expected D factor a^(-2) = %Ps\n",
                        expected_factor);

            I_prime_vect_sigma2 =
                my_H90_vect_2(
                    Labs_cup, Lrel_cup, Lbnr_cup, K, sigma_squared,
                    Ja_vect, p, n);
        }

        for (j=1; j<r_rk+1; ++j) {
            // evaluate cup on j:th basis class (a,J) 
            // DEBUG_PRINT(1, "I'_vect[%d]: %Ps\n\n", j, gel(I_prime_vect, j));
            // DEBUG_PRINT(1, "j=[%d]\n\n", j);
            I_rel = rnfidealabstorel(Lrel_cup, gel(I_prime_vect, j));
            NI = rnfidealnormrel(Lrel_cup, I_rel);
            //DEBUG_PRINT(1, "I, %d to rel\n\n", j);
            
            if (p_int == 3) {
                NIpJ = idealmul(K, gel(gel(Ja_vect, j),2), NI);
            }
            else {
                NIpJ = NI;
            }

            GEN D_H90 = NULL;
            if (diagnostics)
            {
                GEN class_exp = bnfisprincipal0(K, NI, 0);
                D_H90 = my_p_relevant_coordinates(K, class_exp, p);
                pari_printf("\nMASSEY_DIAGNOSTICS norm i=%d, j=%d\n", i, j);
                pari_printf("  matrix row=%d, H90 x-index=%d\n", j, i);
                pari_printf("  I'_absolute = %Ps\n", gel(I_prime_vect, j));
                pari_printf("  norm ideal NI = %Ps\n", NI);
                pari_printf("  norm class D_H90(e_j) = %Ps\n", gtovec(D_H90));
            }

            if (scaling_test)
            {
                GEN I_rel_sigma2 =
                    rnfidealabstorel(
                        Lrel_cup, gel(I_prime_vect_sigma2, j));
                GEN NI_sigma2 =
                    rnfidealnormrel(Lrel_cup, I_rel_sigma2);
                GEN class_exp_sigma2 =
                    bnfisprincipal0(K, NI_sigma2, 0);
                GEN D_sigma2 =
                    my_p_relevant_coordinates(
                        K, class_exp_sigma2, p);
                GEN expected = cgetg(lg(D_H90), t_COL);

                for (long coordinate = 1;
                     coordinate < lg(D_H90); ++coordinate)
                {
                    gel(expected, coordinate) =
                        Fp_mul(expected_factor,
                               gel(D_H90, coordinate), p);
                    if (!equalii(gel(D_sigma2, coordinate),
                                 gel(expected, coordinate)))
                    {
                        pari_printf(
                            "\nMASSEY_DIAGNOSTICS H90 scaling mismatch "
                            "i=%d, j=%d\n", i, j);
                        pari_printf("  sigma = %Ps\n", sigma_cup);
                        pari_printf("  sigma^2 = %Ps\n", sigma_squared);
                        pari_printf("  D_sigma = %Ps\n",
                                    gtovec(D_H90));
                        pari_printf("  D_sigma2 = %Ps\n",
                                    gtovec(D_sigma2));
                        pari_printf("  expected factor = %Ps\n",
                                    expected_factor);
                        pari_printf("  expected = %Ps\n",
                                    gtovec(expected));
                        my_diagnostic_fail(
                            "independent H90 sigma-squared norm class "
                            "does not have the expected scaling");
                    }
                }

                pari_printf(
                    "\nMASSEY_DIAGNOSTICS H90 scaling test i=%d, j=%d\n",
                    i, j);
                pari_printf("  a = 2\n");
                pari_printf("  D_sigma = %Ps\n", gtovec(D_H90));
                pari_printf("  D_sigma2 = %Ps\n", gtovec(D_sigma2));
                pari_printf("  expected = %Ps * D_sigma = %Ps\n",
                            expected_factor, gtovec(expected));
                pari_printf("  OK\n");
            }
            
            
            for (k=1; k<p_rk+1; ++k) {
                DEBUG_PRINT(1, ANSI_COLOR_GREEN "Start: [%d,%d,%d]\n" ANSI_COLOR_RESET, i,j,k);
                // take Artin symbol with resp. to k:th extension
                Labs = gel(gel(K_ext, k), 1);
                Lrel = gel(gel(K_ext, k), 2);
                
                long column = p_rk*(i-1)+k;
                long old_value =
                    my_Artin_symbol(Labs, Lrel, K,
                                    idealred(K,NIpJ), p_int) % p_int;
                gmael2(massey_matrix, j, column) = stoi(old_value);

                if (diagnostics)
                {
                    GEN A_k = gmael(diagnostic_data, k, 2);
                    long new_value = 0;
                    for (long coordinate = 1;
                         coordinate < lg(D_H90); ++coordinate)
                    {
                        long a = umodiu(gel(A_k, coordinate), p_int);
                        long d = umodiu(gel(D_H90, coordinate), p_int);
                        new_value =
                            (new_value + (a * d) % p_int) % p_int;
                    }
                    long direct_value =
                        my_Artin_symbol(Labs, Lrel, K, NI, p_int) % p_int;
                    old_value = (old_value + p_int) % p_int;
                    direct_value = (direct_value + p_int) % p_int;
                    if (old_value != direct_value ||
                        old_value != new_value)
                        my_diagnostic_fail(
                            "Artin symbol does not equal character dot norm class");
                    pari_printf(
                        "  eval k=%d, j=%d: x=%d, y=%d, row=%d, column=%ld, "
                        "class=%Ps, Artin=%ld, dot=%ld, matrix=%Ps, OK\n",
                        k, j, i, k, j, column, gtovec(D_H90),
                        old_value, new_value,
                        gmael2(massey_matrix, j, column));
                }
                DEBUG_PRINT(1, ANSI_COLOR_GREEN "End: [%d,%d,%d]\n\n" ANSI_COLOR_RESET, i,j,k);
                //DEBUG_PRINT(1, "ev_j(x_ix_k): %Ps\n\n", stoi(smodis(my_Artin_symbol(Labs, Lrel, K, NIpJ, p_int, sigma), p_int)));
            }
            
        }
        
    }
    DEBUG_PRINT(0, ANSI_COLOR_MAGENTA "\n-------------------------------------------------------\n %d-fold Massey products of the form < x, x, ..., x, y >\n-------------------------------------------------------\n\n" ANSI_COLOR_RESET, n+1);
    DEBUG_PRINT(1, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, massey_matrix);
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Matrix:  \n\n" ANSI_COLOR_RESET);
    for (j=1; j<nr_row+1; ++j) {
        DEBUG_PRINT(0, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(massey_matrix, j));
    }
    
    
    DEBUG_PRINT(0, ANSI_COLOR_YELLOW "rank: ");
    long mat_rk = FpM_rank((ZM_copy(massey_matrix)), p);
    DEBUG_PRINT(0, ANSI_COLOR_CYAN "%ld\n\n" ANSI_COLOR_RESET, mat_rk);


    if (mat_rk > 0 && n < 3) {
        GEN massey_hnf = FpM_red(hnf((ZM_copy(massey_matrix))),p);
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Hermite normal form:  \n\n" ANSI_COLOR_RESET);
        for (i=1;i<glength(massey_hnf)+1;++i) {
            DEBUG_PRINT(0, ANSI_COLOR_CYAN "%Ps\n\n" ANSI_COLOR_RESET, gel(massey_hnf, i));
        }
        
        DEBUG_PRINT(0, "\n\n");
        //--------------------------------------------------------------------------------
        // Bockstein relations for presentation of Q_2
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Bockstein relations:  \n\n" ANSI_COLOR_RESET);
        char letters[] = "abcdefghijklmnopqr";
        int hnf_r_rk = glength(massey_hnf);
        DEBUG_PRINT(0, "\"\"\"[");
        for (j=1; j<=hnf_r_rk; j++) {
            for (i=1; i<=p_rk; ++i) {
                if (!gequal0(gel(gel(massey_hnf, j), p_rk*(i-1)+i))) {
                    DEBUG_PRINT(0, "%c^%Ps", letters[i-1], gmul(p, gel(gel(massey_hnf, j), p_rk*(i-1)+i)));
                }    
            }
            if (j==hnf_r_rk) {
                DEBUG_PRINT(0, "]\"\"\"\n");
            }
            else {
                DEBUG_PRINT(0, ", ");
            }
        }
        DEBUG_PRINT(0, "\n\n");

        //--------------------------------------------------------------------------------
        // Massey relations for presentation of ZQ_3
        // char letters[] = "abcdefghijklmnopqr";
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Massey relations:  \n\n" ANSI_COLOR_RESET);
        // int hnf_r_rk = glength(massey_hnf);
        DEBUG_PRINT(0, "\"\"\"\"[");
        for (j=1; j<hnf_r_rk+1; j++) {
            for (i=1; i<p_rk+1; ++i) {
                for (k=1; k<p_rk+1; k++) {
                    if (!gequal0(gel(gel(massey_hnf, j), p_rk*(i-1)+k))) {
                        if (p_int==3) {
                            if (i==k) { // there should be no minus sign here since <x,x,x> corresp. to -B(x) (Vogel thesis 1.2.15)
                                DEBUG_PRINT(0, "%c^%Ps", letters[i-1], gmul(p, gel(gel(massey_hnf, j), p_rk*(i-1)+k)));
                            }
                            else {
                                if (i<k) // gives the exponent for [[a_1,a_2], a_1] 
                                {
                                    DEBUG_PRINT(0, "(%c_%c_%c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gneg(gel(gel(massey_hnf, j), p_rk*(i-1)+k)));
                                }
                                else // gives the exponent for [[a_1,a_2], a_2]
                                { 
                                    DEBUG_PRINT(0, "(%c_%c_%c)^%Ps", letters[k-1],letters[i-1], letters[i-1], gel(gel(massey_hnf, j), p_rk*(i-1)+k));
                                }
                            }
                        }   
                        else {
                            if (i<k)
                            {
                                DEBUG_PRINT(0, "(%c_%c_%c)^%Ps", letters[i-1],letters[k-1], letters[i-1], gneg(gel(gel(massey_hnf, j), p_rk*(i-1)+k)));
                            }
                            if (i>k) {
                                DEBUG_PRINT(0, "(%c_%c_%c)^%Ps", letters[k-1],letters[i-1], letters[i-1], gel(gel(massey_hnf, j), p_rk*(i-1)+k));
                            }
                        } 
                    }
                }
            }
            if (j==hnf_r_rk) {
                DEBUG_PRINT(0, "]\"\"\"\"");
            }
            else {
                DEBUG_PRINT(0, ", ");
            }
        }
        DEBUG_PRINT(0, "\n\n");
        //--------------------------------------------------------------------------------
    }
    else {
        DEBUG_PRINT(0, ANSI_COLOR_YELLOW "Massey relations:  \n\n" ANSI_COLOR_RESET);
        DEBUG_PRINT(0, "\"\"\"\"[]\"\"\"\"\n\n");
    }
    

    return mat_rk;
} 

void my_print_massey(GEN K_ext, GEN K, GEN p, int p_int, int p_rk, GEN Ja_vect, int r_rk, GEN best_subgroups) {
    int rk_3_fold, rk_5_fold;
    DEBUG_PRINT(1, "\nSTART: 3-fold\n\n");
    rk_3_fold = my_massey_matrix(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk, best_subgroups, 2);
    DEBUG_PRINT(1, "rk_3_fold: %d\n", rk_3_fold);




    // Higher Massey products
    if (rk_3_fold==0)
    {
        DEBUG_PRINT(1, "\nSTART: 5-fold\n\n");
        rk_5_fold = my_massey_matrix(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk, best_subgroups, 4);
        DEBUG_PRINT(1, "rk_5_fold: %d\n", rk_5_fold);
        if (rk_5_fold==0)
        {
            pari_printf("\nSTART: 7-fold\n\n");
            my_massey_matrix(K_ext, K, p, p_int, p_rk, Ja_vect, r_rk, best_subgroups, 6);
        }
    }
        
    
}
