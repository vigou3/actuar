/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Function to compute the recursive part of the panjer formula
 *  to approximate the aggregate claim amount distribution of
 *  a portfolio over a period.
 *
 *  AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "locale.h"

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))
#define CAD8R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))
#define CAD9R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))))

SEXP panjer(SEXP args)
{
    SEXP p0, p1, fs0, fx0, sfx, a, b, TOL, echo, sfs;
    double *fs, *fx, cumul, constante;
    int r, m, k, x = 1;
	double lastIncr = 1;

    /*  Since we don't know yet what will be the size of the vector fs, 
     *  we'll start by allocating memory for a vector of size 100 with 
     *  S_alloc. That memory space will increase progressively if needed.
     */

    int size = 100;
    fs = (double *) S_alloc(size, sizeof(double));
    for (k = 0; k < size; k++) fs[k] = 0;

    /*  All values received from R are then protected. */

    PROTECT(p0 = coerceVector(CADR(args), REALSXP));
    PROTECT(p1 = coerceVector(CADDR(args), REALSXP));
    PROTECT(fs0 = coerceVector(CADDDR(args), REALSXP));
    PROTECT(fx0 = coerceVector(CAD4R(args), REALSXP));
    PROTECT(sfx = coerceVector(CAD5R(args), REALSXP));
    PROTECT(a = coerceVector(CAD6R(args), REALSXP));
    PROTECT(b = coerceVector(CAD7R(args), REALSXP));
    PROTECT(TOL = coerceVector(CAD8R(args), REALSXP));
    PROTECT(echo = coerceVector(CAD9R(args), LGLSXP));

    fx = REAL(sfx);
    fs[0] = REAL(fs0)[0];
    cumul = REAL(fs0)[0];
    r = length(sfx);
	
	if (LOGICAL(echo)[0]) Rprintf("iteration\t- cumul mpf\t- last increment\n");

    /* (a, b, 0) case (if p0 is NULL) */

    if (isNull(CADR(args)))
    {	
		while (cumul < REAL(TOL)[0] && lastIncr >= R_pow(1-REAL(TOL)[0], 2) )
		{
			if (LOGICAL(echo)[0]) Rprintf("%d\t\t- %.8g\t- %.8g\n", x, cumul, lastIncr);
	    
			if (x >= size) 
			{
			fs = (double *) S_realloc((char *) fs, 2 * size, size, sizeof(double));
			for (k = size; k < 2 * size; k++) fs[k] = 0;
			size = 2 * size;
			}
	    
			m = fmin2(x, r);
	    
			for (k = 1; k <= m; k++)
			{
			fs[x] += ( REAL(a)[0] + REAL(b)[0] * k / x ) * fx[k - 1] * fs[x - k];
			}
	    
			fs[x] = fs[x] / (1 - REAL(a)[0] * REAL(fx0)[0]);
			cumul += fs[x];
			lastIncr = fs[x];
			x++;
		}
		/* raise an error if the Panjer algorithm did not converge*/
		if(cumul < REAL(TOL)[0] && lastIncr < R_pow(1-REAL(TOL)[0], 2) ) 
		error(_("The Panjer algorithm did not converge. You should discretize the claim size distribution on a wider range."));
		
    }
    
    /* (a, b, 1) case (if p0 is non-NULL) */
    
    else
    {
		fx[r] = 0;
		r++;
		constante = (REAL(p1)[0] - (REAL(a)[0] + REAL(b)[0]) * REAL(p0)[0]);
	
		while (cumul < REAL(TOL)[0] && lastIncr >= R_pow(1-REAL(TOL)[0], 2))
		{
			if (LOGICAL(echo)[0]) Rprintf("%d\t\t- %.8g\t- %.8g\n", x, cumul, lastIncr);
	    
			if (x >= size) 
			{
			fs = (double *) S_realloc((char *) fs, 2 * size, size, sizeof(double));
			for (k = size; k < 2 * size; k++) fs[k] = 0;
			size = 2 * size;
			}
	    
			m = fmin2(x, r);
	    
			for (k = 1; k <= m; k++)
			{
			fs[x] += ( REAL(a)[0] + REAL(b)[0] * k / x ) * fx[k - 1] * fs[x - k];
			}
	    
			fs[x] = ( fs[x] + fx[m - 1] * constante ) / (1 - REAL(a)[0] * REAL(fx0)[0]);
			cumul += fs[x];
			lastIncr = fs[x];
			x++;
		}
		/* raise an error if the Panjer algorithm did not converge*/
		if(cumul < REAL(TOL)[0] && lastIncr < R_pow(1-REAL(TOL)[0], 2) ) 
		error(_("The Panjer algorithm did not converge. You should discretize the claim size distribution on a wider range."));
		
    }
    
    /*  Copy of the values to a SEXP which will be returned to R. */
    PROTECT(sfs = allocVector(REALSXP, x));
    memcpy(REAL(sfs), fs, x * sizeof(double));
    
    UNPROTECT(10);
    return(sfs);
}
