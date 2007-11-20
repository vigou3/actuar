/*
 *  Native routines registration, as per "Writing R extensions".
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static const R_ExternalMethodDef ExternalEntries[] = {
    {"do_random", (DL_FUNC) &do_random, -1},
    {"do_dpq", (DL_FUNC) &do_dpq, -1},
    {"do_dpqphtype", (DL_FUNC) &do_dpqphtype, -1},
    {"hierarc", (DL_FUNC) &hierarc, -1},
    {"panjer", (DL_FUNC) &panjer, -1},
    {NULL, NULL, 0}
};

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
}
