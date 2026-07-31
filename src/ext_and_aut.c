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

/**
 * @file ext_and_aut.c
 * @brief Build compatible relative/absolute extension models and automorphisms.
 *
 * The returned package includes relative BNF/BNR data used by later candidate
 * searches.  Its cyclic automorphism is a model generator, not yet the
 * prescribed-character normalization used by the secondary norm.
 */

#include <pari/pari.h>
#include "../headers/ext_and_aut.h"
#include "../headers/misc_functions.h"
#include "../headers/progress.h"

// Debug macros
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define MY_DEBUGLEVEL 0
#define DEBUG_PRINT(level, ...) \
    do { if (MY_DEBUGLEVEL >= (level)) pari_printf(__VA_ARGS__); } while (0)

//------------------------------------------------------------------------------
/*
Generates all extensions of base together with all necessary automorphisms 
*/

GEN my_ext(GEN base, GEN base_clf, GEN p, int p_rk, GEN D_prime_vect) 
{   
    pari_sp av = avma;
    DEBUG_PRINT(1, "\n--------------------------\nStart: my_ext\n--------------------------\n\n");
    GEN x, y, p1, q1, p1red, Lrel, Labs, Lbnr, s_lift_x, cx, sigma, s=pol_x(fetch_user_var("s"));

    x = pol_x(fetch_user_var("x"));
    y = pol_x(fetch_user_var("y"));

    GEN base_ext = cgetg(p_rk+1,t_VEC);

    // DEBUG_PRINT(1, "base l: %ld\n", glength(base_ext));
    // DEBUG_PRINT(1, "Base_clf: %Ps\n\n", base_clf);
    //char filename[100];
    
    int i, j;
    for (i=1; i<p_rk+1; ++i) {
        switch (typ(gel(base_clf, i))) {
            case t_VEC:    p1 = gel(gel(base_clf, i), 1); break;
            case t_POL:    p1 = gel(base_clf, i); break;
        }
        
        q1 = gsubstpol(p1, x, y);
        
        /* Define Lrel/Labs */
        my_progress("reducing the relative defining polynomial "
                    "(rnfpolredbest)");
        p1red = rnfpolredbest(base, mkvec2(q1, D_prime_vect), 0);

        // p1red = q1;
        DEBUG_PRINT(1, "Reduced polynomial for relative extension found\n");
        my_progress("initializing the relative extension (rnfinit)");
        Lrel = rnfinit(base, p1red);

        DEBUG_PRINT(1, "Lrel found\n");
        DEBUG_PRINT(1, "Abs pol: %Ps\n", rnf_get_polabs(Lrel));

        /* The class group of the absolute class field.  At level 2 PARI
         * itself reports the factor base it chose and the timing of its
         * internal phases, which is where a stalled run is diagnosed. */
        my_progress("class group of the degree %ld absolute class field "
                    "(bnfinit)",
                    (long)degpol(rnf_get_polabs(Lrel)));
        Labs = Buchall(rnf_get_polabs(Lrel), nf_FORCE, DEFAULTPREC);
        if (my_progress_level())
        {   /* The cyclic structure is the one number that explains the cost
             * of everything that follows, so it is worth the detour through
             * a temporary string. */
            char *cyc = GENtostr(bnf_get_cyc(Labs));
            my_progress("class group found: cyc = %s", cyc);
            pari_free(cyc);
        }
        my_progress("ray class group of the class field (bnrinit)");
        Lbnr = bnrinit0(Labs,gen_1, 1);

        // May change flag to nf_FORCE
        // Other precisions: MEDDEFAULTPREC, BIGDEFAULTPREC
        //Labs = Buchall_param(rnf_get_polabs(Lrel), 1.5,1.5,4, nf_FORCE, DEFAULTPREC);
        DEBUG_PRINT(1, "Labs found\n");
        DEBUG_PRINT(1, "L_cyc[%d]: %Ps\n", i, bnf_get_cyc(Labs));
        //DEBUG_PRINT(1, "rel_pol[%d]: %Ps\n", i, p1red);
        DEBUG_PRINT(1, "\nabs pol: %Ps\n\n", gsubstpol(rnf_get_polabs(Lrel),y, s));

        s_lift_x = rnfeltup0(Lrel, s, 1);
        //DEBUG_PRINT(1, "\nLift: %Ps\n\n", s_lift_x);
        my_progress("automorphisms of the class field (galoisconj)");
        cx = galoisconj(Labs, NULL);
        //DEBUG_PRINT(1, "\nGalois conj: %Ps\n\n", cx);
        if (glength(cx)==nf_get_degree(bnf_get_nf(Labs)))
        {
            DEBUG_PRINT(1, ANSI_COLOR_GREEN "\n------------------------\nLabs is Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        else {
            DEBUG_PRINT(1, ANSI_COLOR_RED "\n------------------------\nLabs is not Galois over Q\n------------------------\n\n" ANSI_COLOR_RESET);
        }
        

        my_progress("selecting the generator of Gal(L/K) among %ld "
                    "automorphisms", (long)glength(cx));
        for (j = 1; j < glength(cx)+1; ++j)
        {
            if ( (!ZV_equal(algtobasis(Labs,gel(cx, j)), algtobasis(Labs,y))) && ZV_equal(galoisapply(Labs, gel(cx,j), s_lift_x), s_lift_x))
            {
                sigma = algtobasis(Labs, gel(cx, j));
                //DEBUG_PRINT(1, "sigma: %Ps\n\n", sigma);
                break;
            }
        }
        //DEBUG_PRINT(1, ANSI_COLOR_CYAN "---> sigma <--- \n \n" ANSI_COLOR_RESET);
        
    
        gel(base_ext, i) = mkvec4(Labs, Lrel, sigma, Lbnr);

        // sDEBUG_PRINT(1, filename, "/Users/eric/Documents/Matematik/cup-products/large-fields/test/ext_%d", i);
        // writebin(filename, gel(base_ext, i));
        
        
        // my_test_artin_symbol (Labs, Lrel, base, itos(p), sigma);
    }
    DEBUG_PRINT(1, ANSI_COLOR_GREEN "Extensions and generators found\n----------------------------\n\n" ANSI_COLOR_RESET);
    base_ext = gerepilecopy(av, base_ext);
    DEBUG_PRINT(1, "\n--------------------------\nEnd: my_ext\n--------------------------\n\n");
    return base_ext;
}
