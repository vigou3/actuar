#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actuar.h"

static const R_CMethodDef CEntries[]  = {
    {"rpareto", &rpareto, 4, {INTSXP, REALSXP, REALSXP, REALSXP}},
    {NULL, NULL, 0}
}

void R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}

     
	
