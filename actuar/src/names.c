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
    /* DENSITY, CUMULATIVE PROBABILITY AND QUANTILE FUNCTIONS */
    /* One parameter distributions */
    {"dinvexp",		do_dpq1,	1},
    {"pinvexp",		do_dpq1,	2},
    {"qinvexp",		do_dpq1,	3},
    {"minvexp",		do_dpq1,	4},
    {"mexp",		do_dpq1,	5},
    /* Two parameter distributions */
    {"dinvparalogis",	do_dpq2,	1},
    {"pinvparalogis",	do_dpq2,	2},
    {"qinvparalogis",	do_dpq2,	3},
    {"minvparalogis",	do_dpq2,	4},
    {"dinvpareto",	do_dpq2,	5},
    {"pinvpareto",	do_dpq2,	6},
    {"qinvpareto",	do_dpq2,	7},
    {"minvpareto",	do_dpq2,	8},
    {"dinvweibull",	do_dpq2,	9},
    {"pinvweibull",	do_dpq2,	10},
    {"qinvweibull",	do_dpq2,	11},
    {"minvweibull",	do_dpq2,	12},
    {"dllogis",		do_dpq2,	13},
    {"pllogis",		do_dpq2,	14},
    {"qllogis",		do_dpq2,	15},
    {"mllogis",		do_dpq2,	16},
    {"dparalogis",	do_dpq2,	17},
    {"pparalogis",	do_dpq2,	18},
    {"qparalogis",	do_dpq2,	19},
    {"mparalogis",	do_dpq2,	20},
    {"dpareto",		do_dpq2,	21},
    {"ppareto",		do_dpq2,	22},
    {"qpareto",		do_dpq2,	23},
    {"mpareto",		do_dpq2,	24},
    {"dpareto1",	do_dpq2,	25},
    {"ppareto1",	do_dpq2,	26},
    {"qpareto1",	do_dpq2,	27},
    {"mpareto1",	do_dpq2,	28},
    {"dinvgamma",	do_dpq2,	29},
    {"pinvgamma",	do_dpq2,	30},
    {"qinvgamma",	do_dpq2,	31},
    {"minvgamma",	do_dpq2,	32},
    {"dlgamma",	        do_dpq2,	33},
    {"plgamma",	        do_dpq2,	34},
    {"qlgamma",	        do_dpq2,	35},
    {"levinvexp",	do_dpq2,	36}, 
    {"levexp",	        do_dpq2,	37},
    {"mgamma",  	do_dpq2,	38},
    {"mweibull",	do_dpq2,	39},
    {"mlnorm",  	do_dpq2,	40},
    /* Three parameter distributions */
    {"dburr",   	do_dpq3,	1},
    {"pburr",		do_dpq3,	2},
    {"qburr",		do_dpq3,	3},
    {"mburr",		do_dpq3,	4},
    {"dgenpareto",	do_dpq3,	5},
    {"pgenpareto",	do_dpq3,	6},
    {"qgenpareto",	do_dpq3,	7},
    {"mgenpareto",      do_dpq3,        8},
    {"dinvburr",	do_dpq3,	9},
    {"pinvburr",	do_dpq3,	10},
    {"qinvburr",	do_dpq3,	11},
    {"minvburr",	do_dpq3,	12},
    {"dtrgamma",	do_dpq3,	13},
    {"ptrgamma",	do_dpq3,	14},
    {"qtrgamma",	do_dpq3,	15},
    {"mtrgamma",	do_dpq3,	16},
    {"dinvtrgamma",	do_dpq3,	17},
    {"pinvtrgamma",	do_dpq3,	18},
    {"qinvtrgamma",	do_dpq3,	19},
    {"minvtrgamma",	do_dpq3,	20},
    {"levinvparalogis",	do_dpq3,	21},
    {"levinvweibull",	do_dpq3,	22},
    {"levllogis",	do_dpq3,	23},
    {"levparalogis",	do_dpq3,	24},
    {"levpareto",	do_dpq3,	25},
    {"levpareto1",	do_dpq3,	26},
    {"levinvgamma",	do_dpq3,	27},
    {"levgamma",  	do_dpq3,	28},
    {"levweibull",	do_dpq3,	29},
    {"levlnorm",  	do_dpq3,	30},
    /* Four parameter distributions */
    {"dtrbeta",		do_dpq4,	1},
    {"ptrbeta",		do_dpq4,	2},
    {"qtrbeta",		do_dpq4,	3},
    {"mtrbeta",		do_dpq4,	4},
    {"levgenpareto",    do_dpq4,        5},
    {"levburr",		do_dpq4,	6},
    {"levinvburr",	do_dpq4,	7},
    {"levtrgamma",	do_dpq4,	8},
    {"levinvtrgamma",	do_dpq4,	9},
    /* Five parameter distributions */
    {"levtrbeta",	do_dpq5,	1},
	
    /* RANDOM NUMBERS FUNCTIONS */
    /* One parameter distributions */
    {"rinvexp", 	do_random1,	1},
    /* Two parameter distributions */
    {"rinvparalogis", 	do_random2,	1},
    {"rinvpareto", 	do_random2,	2},
    {"rinvweibull",	do_random2,	3},               
    {"rllogis", 	do_random2,	4},
    {"rparalogis", 	do_random2,	5},
    {"rpareto", 	do_random2,	6},
    {"rpareto1", 	do_random2,	7},
    {"rinvgamma", 	do_random2,	8},
    {"rlgamma", 	do_random2,	9},
    {"rinvgauss", 	do_random2,	10},
    /* Three parameter distributions */
    {"rburr", 		do_random3,	1},
    {"rgenpareto", 	do_random3,	2},
    {"rinvburr", 	do_random3,	3},
    {"rtrgamma",	do_random3,	4},
    {"rinvtrgamma",	do_random3,	5},
    /* Four parameter distributions */
    {"rtrbeta", 	do_random4,	1},
    {0, 0, 0}
};
