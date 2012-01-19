/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Table of functions internal to the package. First element is an
 *  argument to one of actuar_do_dpq or actuar_do_random, functions
 *  callable from .External(); second element is the C function
 *  actually called; third element is a code used in the latter.
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
    {"mexp",            actuar_do_dpq1,        1},
    {"dinvexp",         actuar_do_dpq1,        2},
    {"pinvexp",         actuar_do_dpq1,        3},
    {"qinvexp",         actuar_do_dpq1,        4},
    {"minvexp",         actuar_do_dpq1,        5},
    {"mgfexp",          actuar_do_dpq1,        6},
    /* Two parameter distributions */
    {"mgamma",          actuar_do_dpq2,         1},
    {"dinvgamma",       actuar_do_dpq2,         2},
    {"pinvgamma",       actuar_do_dpq2,         3},
    {"qinvgamma",       actuar_do_dpq2,         4},
    {"minvgamma",       actuar_do_dpq2,         5},
    {"dinvparalogis",   actuar_do_dpq2,         6},
    {"pinvparalogis",   actuar_do_dpq2,         7},
    {"qinvparalogis",   actuar_do_dpq2,         8},
    {"minvparalogis",   actuar_do_dpq2,         9},
    {"dinvpareto",      actuar_do_dpq2,        10},
    {"pinvpareto",      actuar_do_dpq2,        11},
    {"qinvpareto",      actuar_do_dpq2,        12},
    {"minvpareto",      actuar_do_dpq2,        13},
    {"dinvweibull",     actuar_do_dpq2,        14},
    {"pinvweibull",     actuar_do_dpq2,        15},
    {"qinvweibull",     actuar_do_dpq2,        16},
    {"minvweibull",     actuar_do_dpq2,        17},
    {"dlgamma",         actuar_do_dpq2,        18},
    {"plgamma",         actuar_do_dpq2,        19},
    {"qlgamma",         actuar_do_dpq2,        20},
    {"mlgamma",         actuar_do_dpq2,        21},
    {"dllogis",         actuar_do_dpq2,        22},
    {"pllogis",         actuar_do_dpq2,        23},
    {"qllogis",         actuar_do_dpq2,        24},
    {"mllogis",         actuar_do_dpq2,        25},
    {"mlnorm",          actuar_do_dpq2,        26},
    {"dparalogis",      actuar_do_dpq2,        27},
    {"pparalogis",      actuar_do_dpq2,        28},
    {"qparalogis",      actuar_do_dpq2,        29},
    {"mparalogis",      actuar_do_dpq2,        30},
    {"dpareto",         actuar_do_dpq2,        31},
    {"ppareto",         actuar_do_dpq2,        32},
    {"qpareto",         actuar_do_dpq2,        33},
    {"mpareto",         actuar_do_dpq2,        34},
    {"dpareto1",        actuar_do_dpq2,        35},
    {"ppareto1",        actuar_do_dpq2,        36},
    {"qpareto1",        actuar_do_dpq2,        37},
    {"mpareto1",        actuar_do_dpq2,        38},
    {"mweibull",        actuar_do_dpq2,        39},
    {"levexp",          actuar_do_dpq2,        40},
    {"levinvexp",       actuar_do_dpq2,        41},
    {"mbeta",           actuar_do_dpq2,        42},
    {"mgfgamma",        actuar_do_dpq2,        43},
    {"mgfnorm",         actuar_do_dpq2,        44},
    {"mgfunif",         actuar_do_dpq2,        45},
    {"mgfinvgamma",     actuar_do_dpq2,        46},
    {"mnorm",           actuar_do_dpq2,        47},
    {"mchisq",          actuar_do_dpq2,        48},
    {"mgfchisq",        actuar_do_dpq2,        49},
    {"minvGauss",       actuar_do_dpq2,        50},
    {"mgfinvGauss",     actuar_do_dpq2,        51},
    {"munif",           actuar_do_dpq2,        52},
    /* Three parameter distributions */
    {"dburr",           actuar_do_dpq3,         1},
    {"pburr",           actuar_do_dpq3,         2},
    {"qburr",           actuar_do_dpq3,         3},
    {"mburr",           actuar_do_dpq3,         4},
    {"dgenpareto",      actuar_do_dpq3,         5},
    {"pgenpareto",      actuar_do_dpq3,         6},
    {"qgenpareto",      actuar_do_dpq3,         7},
    {"mgenpareto",      actuar_do_dpq3,         8},
    {"dinvburr",        actuar_do_dpq3,         9},
    {"pinvburr",        actuar_do_dpq3,        10},
    {"qinvburr",        actuar_do_dpq3,        11},
    {"minvburr",        actuar_do_dpq3,        12},
    {"dinvtrgamma",     actuar_do_dpq3,        13},
    {"pinvtrgamma",     actuar_do_dpq3,        14},
    {"qinvtrgamma",     actuar_do_dpq3,        15},
    {"minvtrgamma",     actuar_do_dpq3,        16},
    {"dtrgamma",        actuar_do_dpq3,        17},
    {"ptrgamma",        actuar_do_dpq3,        18},
    {"qtrgamma",        actuar_do_dpq3,        19},
    {"mtrgamma",        actuar_do_dpq3,        20},
    {"levgamma",        actuar_do_dpq3,        21},
    {"levinvgamma",     actuar_do_dpq3,        22},
    {"levinvparalogis", actuar_do_dpq3,        23},
    {"levinvpareto",    actuar_do_dpq3,        24},
    {"levinvweibull",   actuar_do_dpq3,        25},
    {"levlgamma",       actuar_do_dpq3,        26},
    {"levllogis",       actuar_do_dpq3,        27},
    {"levlnorm",        actuar_do_dpq3,        28},
    {"levparalogis",    actuar_do_dpq3,        29},
    {"levpareto",       actuar_do_dpq3,        30},
    {"levpareto1",      actuar_do_dpq3,        31},
    {"levweibull",      actuar_do_dpq3,        32},
    {"levbeta",         actuar_do_dpq3,        33},
    {"levchisq",        actuar_do_dpq3,        34},
    {"levinvGauss",     actuar_do_dpq3,        35},
    {"levunif",         actuar_do_dpq3,        36},
    /* Four parameter distributions */
    {"dtrbeta",         actuar_do_dpq4,        1},
    {"ptrbeta",         actuar_do_dpq4,        2},
    {"qtrbeta",         actuar_do_dpq4,        3},
    {"mtrbeta",         actuar_do_dpq4,        4},
    {"levburr",         actuar_do_dpq4,        5},
    {"levgenpareto",    actuar_do_dpq4,        6},
    {"levinvburr",      actuar_do_dpq4,        7},
    {"levinvtrgamma",   actuar_do_dpq4,        8},
    {"levtrgamma",      actuar_do_dpq4,        9},
    {"dgenbeta",        actuar_do_dpq4,        10},
    {"pgenbeta",        actuar_do_dpq4,        11},
    {"qgenbeta",        actuar_do_dpq4,        12},
    {"mgenbeta",        actuar_do_dpq4,        13},
    /* Five parameter distributions */
    {"levtrbeta",       actuar_do_dpq5,        1},
    {"levgenbeta",      actuar_do_dpq5,        2},
    /* Phase-type distributions */
    {"dphtype",         actuar_do_dpqphtype2,  1},
    {"pphtype",         actuar_do_dpqphtype2,  2},
    {"mphtype",         actuar_do_dpqphtype2,  3},
    {"mgfphtype",       actuar_do_dpqphtype2,  4},

    /* RANDOM NUMBERS FUNCTIONS */
    /* One parameter distributions */
    {"rinvexp",         actuar_do_random1,     1},
    /* Two parameter distributions */
    {"rinvgamma",       actuar_do_random2,     1},
    {"rinvparalogis",   actuar_do_random2,     2},
    {"rinvpareto",      actuar_do_random2,     3},
    {"rinvweibull",     actuar_do_random2,     4},
    {"rlgamma",         actuar_do_random2,     5},
    {"rllogis",         actuar_do_random2,     6},
    {"rparalogis",      actuar_do_random2,     7},
    {"rpareto",         actuar_do_random2,     8},
    {"rpareto1",        actuar_do_random2,     9},
    /* Three parameter distributions */
    {"rburr",           actuar_do_random3,     1},
    {"rgenpareto",      actuar_do_random3,     2},
    {"rinvburr",        actuar_do_random3,     3},
    {"rinvtrgamma",     actuar_do_random3,     4},
    {"rtrgamma",        actuar_do_random3,     5},
    /* Four parameter distributions */
    {"rtrbeta",         actuar_do_random4,     1},
    {"rgenbeta",        actuar_do_random4,     2},
    /* Phase-type distributions */
    {"rphtype",         actuar_do_randomphtype2, 1},
    {0, 0, 0}
};
