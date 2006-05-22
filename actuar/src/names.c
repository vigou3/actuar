/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Table of functions internal to the package. First element is an
 *  argument to one of do_math or do_random, functions callable from
 *  .External(); second element is the C function actually called;
 *  third element is a code used in the latter.
 *
 *  Idea taken from R sources (see .../src/main/names.c).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <Rinternals.h>
#include "actuar.h"

FUNTAB fun_tab[] = { 
    {"rbetatrans", do_random4, 0},
    {"rparetogen", do_random3, 0},
    {"rburr", do_random3, 1},
    {"rinverseburr", do_random3, 2},
    {"rloglogistic", do_random2, 0},
    {"rparalogistic", do_random2, 1},
    {"rpareto", do_random2, 2},
    {"rinverseparalogistic", do_random2, 3},
    {"rinversepareto", do_random2, 4},
    {"rinverseweibull", do_random2, 5},
    {"rspareto", do_random2, 6},
    {"rinverseexp", do_random1, 0},
    {0, 0, 0}
};
