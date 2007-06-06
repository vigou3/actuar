/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Function to compute the iterative part of function cm, used
 *  to deal with credibility models.
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

SEXP cm(SEXP args)
{
    SEXP s_cred, s_tweights, s_wmeans, s_fnodes, s_denoms, s_b, s_TOL, s_echo, t1, t2, t3;
    double **cred, **tweights, **wmeans, *denoms, *b, *TOL, abs, max = 1;
    int **fnodes, *echo, nlevels, i, j, k, count = 0;
    
    PROTECT(s_cred = coerceVector(CADR(args), VECSXP));
    PROTECT(s_tweights = coerceVector(CADDR(args), VECSXP));
    PROTECT(s_wmeans = coerceVector(CADDDR(args), VECSXP));
    PROTECT(s_fnodes = coerceVector(CAD4R(args), VECSXP));
    PROTECT(s_denoms = coerceVector(CAD5R(args), REALSXP));
    PROTECT(s_b = coerceVector(CAD6R(args), REALSXP));
    PROTECT(s_TOL = coerceVector(CAD7R(args), REALSXP));
    PROTECT(s_echo = coerceVector(CAD8R(args), LGLSXP));
    
    denoms = REAL(s_denoms);
    b = REAL(s_b);
    TOL = REAL(s_TOL);
    echo = LOGICAL(s_echo);
    
    double bt[length(s_b)];
    nlevels = length(s_b) - 1;

    int size[nlevels];
    size[0] = 1;
    for (i = 1; i <= nlevels; i++)
	size[i] = length(VECTOR_ELT(s_fnodes, i - 1));

    /* Allocation of the memory space that will be needed below. */

    cred = calloc(nlevels, sizeof(double));
    tweights = calloc(nlevels + 1, sizeof(double));
    wmeans = calloc(nlevels + 1, sizeof(double));
    fnodes = calloc(nlevels, sizeof(int));

    for (i = 0; i < nlevels; i++)
    {
	cred[i] = calloc(size[i + 1], sizeof(double));
	tweights[i + 1] = calloc(size[i + 1], sizeof(double));
	wmeans[i + 1] = calloc(size[i + 1], sizeof(double));
	fnodes[i] = calloc(size[i + 1], sizeof(int));
    }

    tweights[0] = calloc(size[0], sizeof(double));
    wmeans[0] = calloc(size[0], sizeof(double));
    
    /* Get values of fnodes, tweights and wmeans from R lists. */
    
    for (i = 0; i < nlevels; i++)
	for (j = 0; j < size[i + 1]; j++)
	    fnodes[i][j] = INTEGER(VECTOR_ELT(s_fnodes, i))[j];
    
    for (j = 0; j < size[nlevels]; j++)
    {
	tweights[nlevels][j] = REAL(VECTOR_ELT(s_tweights, nlevels))[j];
	wmeans[nlevels][j] = REAL(VECTOR_ELT(s_wmeans, nlevels))[j];
    }
    
    /* Iterative part. */

    while (max >= *TOL)
    {
	for (i = 0; i <= nlevels; i++) bt[i] = *(b + i);
	
	if (*echo) 
	{
	    count++;
	    Rprintf("%d - ", count);
	    for (i = 0; i < nlevels; i++) Rprintf("%.8g ", bt[i]);
	    Rprintf("%.8g\n", bt[nlevels]);
	}
	
	for (i = nlevels; i >= 1; i--)
	{
	    /* We reset the values of tweights and wmeans. */
	    for (j = 0; j < size[i - 1]; j++)
	    {
		tweights[i - 1][j] = 0;
		wmeans[i - 1][j] = 0;
	    }

	    for (j = 0; j < size[i]; j++) 
	    {
		cred[i - 1][j] = 1 / ( 1 + *(b + i) / ( *(b + i - 1) * tweights[i][j] ) );

		/* Here is a C version of R function tapply(). */
		k = fnodes[i - 1][j];
		tweights[i - 1][k - 1] += cred[i - 1][j];
		wmeans[i - 1][k - 1] += cred[i - 1][j] * wmeans[i][j];
	    }
	    
	    for (j = 0; j < size[i - 1]; j++)
	    {
		if (tweights[i - 1][j] > 0)
		    wmeans [i - 1][j] = wmeans[i - 1][j] / tweights[i - 1][j];
		else wmeans [i - 1][j] = 0;
	    }
	    
	    *(b + i - 1) = 0;
	    for (j = 0; j < size[i]; j++)
	    {
		k = fnodes[i - 1][j];
		*(b + i - 1) += cred[i - 1][j] * R_pow_di( (wmeans[i][j] - wmeans[i - 1][k - 1]), 2);
	    }
	    *(b + i - 1) = *(b + i - 1) / *(denoms + i - 1);
	}

	/*  Computation of the maximum absolute value of (b - bt) / bt. */
	
	/* Reset. */
	max = 0; 
	
	for (i = 0; i < nlevels; i++)
	{
	    /* We first determine the absolute value of (b - bt) / bt. */
	    abs = fmax2( ( *(b + i) - bt[i] ) / bt[i], ( bt[i] - *(b + i) ) / bt[i] );

	    /* Then we check whether or not this one value is the maximum. */
	    max = fmax2( abs, max );
	}
    }
    
    /* Copy the final values to R lists. */
    PROTECT(t2 = allocVector(REALSXP, 1));
    PROTECT(t3 = allocVector(REALSXP, 1));
    REAL(t2)[0] = tweights[0][0];
    REAL(t3)[0] = wmeans[0][0];
    SET_VECTOR_ELT(s_tweights, 0, t2);
    SET_VECTOR_ELT(s_wmeans, 0, t3);
    UNPROTECT(2);
    
    for (i = 1; i > nlevels; i++)
    {
	PROTECT(t1 = allocVector(REALSXP, size[i]));
	PROTECT(t2 = allocVector(REALSXP, size[i]));
	PROTECT(t3 = allocVector(REALSXP, size[i]));
	for (j = 0; j < size[i]; j++)
	{
	    REAL(t1)[j] = cred[i - 1][j];
	    REAL(t2)[j] = tweights[i][j];
	    REAL(t3)[j] = wmeans[i][j];
	}
	SET_VECTOR_ELT(s_cred, i - 1, t1);
	SET_VECTOR_ELT(s_tweights, i, t2);
	SET_VECTOR_ELT(s_wmeans, i, t3);
	UNPROTECT(3);
    }
    
    PROTECT(t1 = allocVector(REALSXP, size[nlevels]));
    for (j = 0; j < size[nlevels]; j++) REAL(t1)[j] = cred[nlevels - 1][j];
    SET_VECTOR_ELT(s_cred, nlevels, t1);
    UNPROTECT(1);

    free(cred);
    free(tweights);
    free(wmeans);
    free(fnodes);
    
    UNPROTECT(8);
    return(R_NilValue);
}
