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
    /* DENSITY, CUMULATIVE PROBABILITY AND QUANTILE FUNCTIONS,
     * RAW AND LIMITED MOMENTS */
    /* One parameter distributions */
    {"mexp",		do_dpq1,	1},
    {"dinvexp",		do_dpq1,	2},
    {"pinvexp",		do_dpq1,	3},
    {"qinvexp",		do_dpq1,	4},
    {"minvexp",		do_dpq1,	5},
    /* Two parameter distributions */
    {"mgamma",  	do_dpq2,	 1},
    {"dinvgamma",	do_dpq2,	 2},
    {"pinvgamma",	do_dpq2,	 3},
    {"qinvgamma",	do_dpq2,	 4},
    {"minvgamma",	do_dpq2,	 5},
    {"dinvparalogis",	do_dpq2,	 6},
    {"pinvparalogis",	do_dpq2,	 7},
    {"qinvparalogis",	do_dpq2,	 8},
    {"minvparalogis",	do_dpq2,	 9},
    {"dinvpareto",	do_dpq2,	10},
    {"pinvpareto",	do_dpq2,	11},
    {"qinvpareto",	do_dpq2,	12},
    {"minvpareto",	do_dpq2,	13},
    {"dinvweibull",	do_dpq2,	14},
    {"pinvweibull",	do_dpq2,	15},
    {"qinvweibull",	do_dpq2,	16},
    {"minvweibull",	do_dpq2,	17},
    {"dlgamma",	        do_dpq2,	18},
    {"plgamma",	        do_dpq2,	19},
    {"qlgamma",	        do_dpq2,	20},
    {"mlgamma",	        do_dpq2,	21},
    {"dllogis",		do_dpq2,	22},
    {"pllogis",		do_dpq2,	23},
    {"qllogis",		do_dpq2,	24},
    {"mllogis",		do_dpq2,	25},
    {"mlnorm",  	do_dpq2,	26},
    {"dparalogis",	do_dpq2,	27},
    {"pparalogis",	do_dpq2,	28},
    {"qparalogis",	do_dpq2,	29},
    {"mparalogis",	do_dpq2,	30},
    {"dpareto",		do_dpq2,	31},
    {"ppareto",		do_dpq2,	32},
    {"qpareto",		do_dpq2,	33},
    {"mpareto",		do_dpq2,	34},
    {"dpareto1",	do_dpq2,	35},
    {"ppareto1",	do_dpq2,	36},
    {"qpareto1",	do_dpq2,	37},
    {"mpareto1",	do_dpq2,	38},
    {"mweibull",	do_dpq2,	39},
    {"levexp",	        do_dpq2,	40},
    {"levinvexp",	do_dpq2,	41},
    {"mbeta",	        do_dpq2,	42},
    /* Three parameter distributions */
    {"dburr",   	do_dpq3,	 1},
    {"pburr",		do_dpq3,	 2},
    {"qburr",		do_dpq3,	 3},
    {"mburr",		do_dpq3,	 4},
    {"dgenpareto",	do_dpq3,	 5},
    {"pgenpareto",	do_dpq3,	 6},
    {"qgenpareto",	do_dpq3,	 7},
    {"mgenpareto",      do_dpq3,         8},
    {"dinvburr",	do_dpq3,	 9},
    {"pinvburr",	do_dpq3,	10},
    {"qinvburr",	do_dpq3,	11},
    {"minvburr",	do_dpq3,	12},
    {"dinvtrgamma",	do_dpq3,	13},
    {"pinvtrgamma",	do_dpq3,	14},
    {"qinvtrgamma",	do_dpq3,	15},
    {"minvtrgamma",	do_dpq3,	16},
    {"dtrgamma",	do_dpq3,	17},
    {"ptrgamma",	do_dpq3,	18},
    {"qtrgamma",	do_dpq3,	19},
    {"mtrgamma",	do_dpq3,	20},
    {"levgamma",  	do_dpq3,	21},
    {"levinvgamma",	do_dpq3,	22},
    {"levinvparalogis",	do_dpq3,	23},
    {"levinvpareto",	do_dpq3,	24},
    {"levinvweibull",	do_dpq3,	25},
    {"levlgamma",  	do_dpq3,	26},
    {"levllogis",	do_dpq3,	27},
    {"levlnorm",  	do_dpq3,	28},
    {"levparalogis",	do_dpq3,	29},
    {"levpareto",	do_dpq3,	30},
    {"levpareto1",	do_dpq3,	31},
    {"levweibull",	do_dpq3,	32},
    {"levbeta",	        do_dpq3,	33},
	{"mgfgamma",    do_dpq3,    34},
    /* Four parameter distributions */
    {"dtrbeta",		do_dpq4,	1},
    {"ptrbeta",		do_dpq4,	2},
    {"qtrbeta",		do_dpq4,	3},
    {"mtrbeta",		do_dpq4,	4},
    {"levburr",		do_dpq4,	5},
    {"levgenpareto",    do_dpq4,        6},
    {"levinvburr",	do_dpq4,	7},
    {"levinvtrgamma",	do_dpq4,	8},
    {"levtrgamma",	do_dpq4,	9},
    {"dgenbeta",	do_dpq4,        10},
    {"pgenbeta",	do_dpq4,        11},
    {"qgenbeta",	do_dpq4,        12},
    {"mgenbeta",	do_dpq4,        13},
    /* Five parameter distributions */
    {"levtrbeta",	do_dpq5,	1},
    {"levgenbeta",	do_dpq5,	2},

    /* RANDOM NUMBERS FUNCTIONS */
    /* One parameter distributions */
    {"rinvexp", 	do_random1,	1},
    /* Two parameter distributions */
    {"rinvgamma", 	do_random2,	1},
    {"rinvparalogis", 	do_random2,	2},
    {"rinvpareto", 	do_random2,	3},
    {"rinvweibull",	do_random2,	4},
    {"rlgamma", 	do_random2,	5},
    {"rllogis", 	do_random2,	6},
    {"rparalogis", 	do_random2,	7},
    {"rpareto", 	do_random2,	8},
    {"rpareto1", 	do_random2,	9},
    /* Three parameter distributions */
    {"rburr", 		do_random3,	1},
    {"rgenpareto", 	do_random3,	2},
    {"rinvburr", 	do_random3,	3},
    {"rinvtrgamma",	do_random3,	4},
    {"rtrgamma",	do_random3,	5},
    /* Four parameter distributions */
    {"rtrbeta", 	do_random4,	1},
    {"rgenbeta", 	do_random4,	2},
    {0, 0, 0}
};
