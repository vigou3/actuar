/*
 *  Native routines registration, as per "Writing R extensions" and
 *  definition of native interface to one routine exported by package
 *  expint.
 *
 */

#include <Rinternals.h>
#include <R_ext/Visibility.h>
#include <R_ext/Rdynload.h>
#include <expintAPI.h>		/* optional; for expressiveness */
#include "actuar.h"

static const R_ExternalMethodDef ExternalEntries[] = {
    {"actuar_do_random", (DL_FUNC) &actuar_do_random, -1},
    {"actuar_do_randomphtype", (DL_FUNC) &actuar_do_randomphtype, -1},
    {"actuar_do_dpq", (DL_FUNC) &actuar_do_dpq, -1},
    {"actuar_do_dpqphtype", (DL_FUNC) &actuar_do_dpqphtype, -1},
    {"actuar_do_hierarc", (DL_FUNC) &actuar_do_hierarc, -1},
    {"actuar_do_panjer", (DL_FUNC) &actuar_do_panjer, -1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);

    /* native interface to routine from package expint */
    actuar_gamma_inc = (double(*)(double,double)) R_GetCCallable("expint", "gamma_inc");
}
