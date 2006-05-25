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
    {"ppareto",			do_dpq2,	2},
    {"qpareto",			do_dpq2,	3},
    {"dpareto1",		do_dpq2,	4},
    {"ppareto1",		do_dpq2,	5},
    {"qpareto1",		do_dpq2,	6},

    {"dburr",	        	do_dpq3,	0},
    {"pburr",		        do_dpq3,	1},
    {"qburr",		        do_dpq3,	2},
	
    {"dtrbeta",	        	do_dpq4,	0},
    {"ptrbeta",		        do_dpq4,	1},
    {"qtrbeta",		        do_dpq4,	2},

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
