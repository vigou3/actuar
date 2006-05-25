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
    {"dpareto",			do_dpq2,	1},
    {"qpareto",			do_dpq2,	2},
    {"rpareto",			do_dpq2,	3},
	
    {"rexp", 		        do_random1,	0},

    {"rllogis", 		do_random2,	0},
    {"rparalogis", 		do_random2,	1},
    {"rpareto", 		do_random2,	2},
    {"riparalogis",     	do_random2,	3},
    {"ripareto", 		do_random2,	4},
    {"rlgompertz",      	do_random2,	5},
    {"rpareto1", 		do_random2,	6},

    {"rgenpareto", 		do_random3,	0},
    {"rburr", 			do_random3,	1},
    {"riburr", 	        	do_random3,	2},

    {"rtrbeta", 		do_random4,	0},
    {0, 0, 0}
};
