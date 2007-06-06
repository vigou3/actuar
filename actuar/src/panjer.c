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

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))
#define CAD8R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))
#define CAD9R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))))

SEXP panjer(SEXP args)
{
    SEXP sp0, sp1, sfs0, sfx0, sfx, sa, sb, sTOL, secho, sfs;
    double *fs, *p0, *p1, *fs0, *fx0, *fx, *a, *b, *TOL, *fs2, cumul, constante;
    int *echo, r, m, k, x = 1;

    /*  Calloc allows us to allocate memory for the vector fs, which 
     *  will be used to stock the values of the claim amount in the loop 
     *  below. Since we don't know how many iterations will be needed, we 
     *  first start by setting the size to 100, size that will increase 
     *  progressively if needed.
     */

    int size = 100;
    fs = Calloc(size, double);
    for (k = 0; k < size; k++) fs[k] = 0;

    /*  All values received from R are then associated with a numeric
     *  or logical pointer.
     */

    PROTECT(sp0 = coerceVector(CADR(args), REALSXP));
    PROTECT(sp1 = coerceVector(CADDR(args), REALSXP));
    PROTECT(sfs0 = coerceVector(CADDDR(args), REALSXP));
    PROTECT(sfx0 = coerceVector(CAD4R(args), REALSXP));
    PROTECT(sfx = coerceVector(CAD5R(args), REALSXP));
    PROTECT(sa = coerceVector(CAD6R(args), REALSXP));
    PROTECT(sb = coerceVector(CAD7R(args), REALSXP));
    PROTECT(sTOL = coerceVector(CAD8R(args), REALSXP));
    PROTECT(secho = coerceVector(CAD9R(args), LGLSXP));

    p0 = REAL(sp0);    
    p1 = REAL(sp1);
    fs0 = REAL(sfs0);
    fx0 = REAL(sfx0);
    fx = REAL(sfx);
    a = REAL(sa);
    b = REAL(sb);
    TOL = REAL(sTOL);
    echo = LOGICAL(secho);

    fs[0] = *fs0;
    cumul = *fs0;
    r = length(sfx);


    /* (a, b, 0) case (if p0 is NULL) */

    if (isNull(CADR(args)))
    {	
	while (cumul < *TOL)
	{
	    if (*echo) Rprintf("%d - %.8g\n", x, cumul);
	    
	    if (x >= size) 
	    {
		size = 2 * size;
		fs = Realloc(fs, size, double);
		for (k = size / 2; k < size; k++) fs[k] = 0;
	    }
	    
	    m = fmin2(x, r);
	    
	    for (k = 1; k <= m; k++)
	    {
		fs[x] += ( *a + *b * k / x ) * *(fx + k - 1) * fs[x - k];
	    }
	    
	    fs[x] = fs[x] / (1 - *a * *fx0);
	    cumul += fs[x];
	    x++;
	}
    }
    
    
    /* (a, b, 1) case (if p0 is non-NULL) */
    
    else
    {
	*(fx + r) = 0;
	r++;
	constante = (*p1 - (*a + *b) * *p0);
	
	while (cumul < *TOL)
	{
	    if (*echo) Rprintf("%d - %.8g\n", x, cumul);
	    
	    if (x >= size) 
	    {
		size = 2 * size;
		fs = Realloc(fs, size, double);
		for (k = size / 2; k < size; k++) fs[k] = 0;
	    }
	    
	    m = fmin2(x, r);
	    
	    for (k = 1; k <= m; k++)
	    {
		fs[x] += ( *a + *b * k / x ) * *(fx + k - 1) * fs[x - k];
	    }
	    
	    fs[x] = ( fs[x] + *(fx + m - 1) * constante ) / (1 - *a * *fx0);
	    cumul += fs[x];
	    x++;
	}
    }
    
    /*  A new variable of the correct length, fs2, is created to get rid of
     *  the zeros and the extra space that was accorded to fs. That new 
     *  variable points towards a SEXP that will be returned to R.
     */

    PROTECT(sfs = allocVector(REALSXP,x));
    fs2 = REAL(sfs);
    
    for (k = 0; k < x; k++) *(fs2 + k) = fs[k];
    
    Free(fs);

    UNPROTECT(10);
    return(sfs);
}
