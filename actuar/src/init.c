#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static R_NativePrimitiveArgType rpareto_t[4] = {INTSXP, REALSXP, REALSXP, REALSXP};

static const R_CMethodDef CEntries[] = {
    {"rpareto", (DL_FUNC) &rpareto, 4, rpareto_t},
    {NULL, NULL, 0}
};

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
