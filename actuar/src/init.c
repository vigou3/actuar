#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static R_NativePrimitiveArgType rpareto_t[6] = {INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};

static const R_CMethodDef CEntries[] = {
    {"rpareto_sym", (DL_FUNC) &rpareto, 6, rpareto_t},
    {NULL, NULL, 0}
};

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}
