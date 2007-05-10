/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Function to compute the recursive part of the panjer formula
 *  to approximate the aggregate claim amount distribution of
 *  a portfolio over a period.
 *
 *  AUTHOR: Tommy Ouellet
 */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "locale.h"
#include "dpq.h"

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))
#define CAD8R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))
#define CAD9R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))))

SEXP panjer(SEXP args)
{
    SEXP p0_1, p1_1, fs0_1, fx0_1, fx_1, a_1, b_1, TOL_1, echo_1, fs_1;
    double *fs, *p0_2, *p1_2, *fs0_2, *fx0_2, *fx_2, *a_2, *b_2, *TOL_2, *fs_2, cumul, constante;
    int *echo_2, r, m, k, x = 1, size=100;
    
    fs = malloc(size * sizeof(double));
    for (k = 0; k < size; k++) fs[k] = 0;

    PROTECT(p0_1 = AS_NUMERIC(CADR(args)));
    PROTECT(p1_1 = AS_NUMERIC(CADDR(args)));
    PROTECT(fs0_1 = AS_NUMERIC(CADDDR(args)));
    PROTECT(fx0_1 = AS_NUMERIC(CAD4R(args)));
    PROTECT(fx_1 = AS_NUMERIC(CAD5R(args)));
    PROTECT(a_1 = AS_NUMERIC(CAD6R(args)));
    PROTECT(b_1 = AS_NUMERIC(CAD7R(args)));
    PROTECT(TOL_1 = AS_NUMERIC(CAD8R(args)));
    PROTECT(echo_1 = AS_LOGICAL(CAD9R(args)));

    p0_2 = NUMERIC_POINTER(p0_1);
    p1_2 = NUMERIC_POINTER(p1_1);
    fs0_2 = NUMERIC_POINTER(fs0_1);
    fx0_2 = NUMERIC_POINTER(fx0_1);
    fx_2 = NUMERIC_POINTER(fx_1);
    a_2 = NUMERIC_POINTER(a_1);
    b_2 = NUMERIC_POINTER(b_1);
    TOL_2 = NUMERIC_POINTER(TOL_1);
    echo_2 = LOGICAL_POINTER(echo_1);

    fs[0] = *fs0_2;
    cumul = *fs0_2;
    r = length(fx_1);


    /* Classe (a,b,0) */

    if (*p0_2 == -1)
    {	
	while (*TOL_2 > cumul)
	{
	    if (*echo_2) Rprintf("%d - %.8f\n", x, cumul);
	    
	    if (x >= size) 
	    {
		size = 2 * size;
		fs = realloc(fs, size * sizeof(double));
		for (k = size/2; k < size; k++) fs[k] = 0;
	    }
	    
	    if (x > r) m = r; else m = x;
	    
	    for (k = 1; k <= m; k++)
	    {
		fs[x] = fs[x] + ( *a_2 + *b_2 * k / x ) * *(fx_2+k-1) * fs[x-k];
	    }
	    
	    fs[x] = fs[x] / (1 - *a_2 * *fx0_2);
	    cumul = cumul + fs[x];
	    x++;
	}
    }
    
    
    /* Classe (a,b,1) */
    
    else
    {
	*(fx_2+r) = 0;
	r++;
	constante = (*p1_2 - (*a_2 + *b_2) * *p0_2);
	
	while (*TOL_2 > cumul)
	{
	    if (*echo_2) Rprintf("%d - %.8f\n", x, cumul);
	    
	    if (x >= size) 
	    {
		size = 2 * size;
		fs = realloc(fs, size * sizeof(double));
		for (k = size/2; k < size; k++) fs[k] = 0;
	    }
	    
	    if (x > r) m= r; else m = x;
	    
	    for (k = 1; k <= m; k++)
	    {
		fs[x] = fs[x] + ( *a_2 + *b_2 * k / x ) * *(fx_2+k-1) * fs[x-k];
	    }
	    
	    fs[x] = ( fs[x] + *(fx_2+m-1) * constante ) / (1 - *a_2 * *fx0_2);
	    cumul = cumul + fs[x];
	    x++;
	}
    }
    
    PROTECT(fs_1 = NEW_NUMERIC(size));
    fs_2 = NUMERIC_POINTER(fs_1);
    
    for (k = 0; k < size; k++) *(fs_2+k) = fs[k];
    
    free(fs);
    fs = NULL;

    UNPROTECT(10);
    return(fs_1);
}

