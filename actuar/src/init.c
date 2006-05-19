/*
 *  Native routines registration, as per "Writing R extensions".
 */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static const R_ExternalMethodDef ExternalEntries[] = {
    {"do_random", (DL_FUNC) &do_random, -1},
    {NULL, NULL, 0}
};

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
}
