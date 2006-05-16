#include <R.h>
#include <Rinternals.h>

const R_CMethodDef CEntries[]  = {
    {"rpareto", &rpareto, 4, {INTSXP, REALSXP, REALSXP, REALSXP}},
    {NULL, NULL, 0}
}

R_init_actuar(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
}

     
	
