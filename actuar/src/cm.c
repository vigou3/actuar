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

#define abs(x) (x >= 0 ? x : -x)
#define tozero(x) (x > R_pow_di(REAL(TOL)[0], 2) ? x : 0)
#define weights(i, j) ( cred[i - 1][j] != 0 ? cred[i - 1][j] : tweights[i][j] )

SEXP toSEXP(double *x, int size)
{
    SEXP ans = allocVector(REALSXP, size);
    memcpy(REAL(ans), x, size * sizeof(double));
    return ans;
}

SEXP cm(SEXP args)
{
    SEXP s_cred, s_tweights, s_wmeans, s_fnodes, denoms, b, TOL, echo;
    double **cred, **tweights, **wmeans, max = 1;
    int **fnodes, nlevels, i, j, k, count = 0;
    
    PROTECT(s_cred = coerceVector(CADR(args), VECSXP));
    PROTECT(s_tweights = coerceVector(CADDR(args), VECSXP));
    PROTECT(s_wmeans = coerceVector(CADDDR(args), VECSXP));
    PROTECT(s_fnodes = coerceVector(CAD4R(args), VECSXP));
    PROTECT(denoms = coerceVector(CAD5R(args), REALSXP));
    PROTECT(b = coerceVector(CAD6R(args), REALSXP));
    PROTECT(TOL = coerceVector(CAD7R(args), REALSXP));
    PROTECT(echo = coerceVector(CAD8R(args), LGLSXP));
    
    double bt[length(b)];
    nlevels = length(b) - 1;

    int size[nlevels]; size[0] = 1;
    for (i = 1; i <= nlevels; i++)
	size[i] = length(VECTOR_ELT(s_fnodes, i - 1));

    /* Allocation of the memory space that will be needed below. */

    cred = (double **) S_alloc(nlevels, sizeof(double));
    tweights = (double **) S_alloc(nlevels + 1, sizeof(double));
    wmeans = (double **) S_alloc(nlevels + 1, sizeof(double));
    fnodes = (int **) S_alloc(nlevels, sizeof(int));

    for (i = 0; i < nlevels; i++)
    {
	cred[i] = (double *) S_alloc(size[i + 1], sizeof(double));
	tweights[i + 1] = (double *) S_alloc(size[i + 1], sizeof(double));
	wmeans[i + 1] = (double *) S_alloc(size[i + 1], sizeof(double));
	fnodes[i] = (int *) S_alloc(size[i + 1], sizeof(int));
    }

    tweights[0] = (double *) S_alloc(size[0], sizeof(double));
    wmeans[0] = (double *) S_alloc(size[0], sizeof(double));
    
    /* Get values of fnodes, tweights and wmeans from R lists. */
    
    for (i = 0; i < nlevels; i++)
	memcpy(fnodes[i], INTEGER(VECTOR_ELT(s_fnodes, i)), size[i + 1] * sizeof(int));
    memcpy(tweights[nlevels], REAL(VECTOR_ELT(s_tweights, nlevels)), size[nlevels] * sizeof(double));
    memcpy(wmeans[nlevels], REAL(VECTOR_ELT(s_wmeans, nlevels)), size[nlevels] * sizeof(double));

    /* Iterative part. */

    while (max >= REAL(TOL)[0])
    {
	for (i = 0; i <= nlevels; i++) bt[i] = REAL(b)[i];
	
	if (LOGICAL(echo)[0]) 
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
		cred[i - 1][j] = 1 / ( 1 + REAL(b)[i] / ( REAL(b)[i - 1] * tweights[i][j] ) );

		/* Here is a C version of R function tapply(). */
		k = fnodes[i - 1][j];
		tweights[i - 1][k - 1] += weights(i, j);
		wmeans[i - 1][k - 1] += weights(i, j) * wmeans[i][j];
	    }
	    
	    for (j = 0; j < size[i - 1]; j++)
	    {
		if (tweights[i - 1][j] > 0)
		    wmeans [i - 1][j] = wmeans[i - 1][j] / tweights[i - 1][j];
		else wmeans [i - 1][j] = 0;
	    }
	    
	    if (REAL(b)[i - 1])
	    {
		REAL(b)[i - 1] = 0;
		for (j = 0; j < size[i]; j++)
		{
		    k = fnodes[i - 1][j];
		    REAL(b)[i - 1] += weights(i, j) * R_pow_di( (wmeans[i][j] - wmeans[i - 1][k - 1]), 2);
		}
		REAL(b)[i - 1] = REAL(b)[i - 1] / REAL(denoms)[i - 1];
	    }
	}
	
	/*  Computation of the maximum absolute value of (b - bt) /
	 *  bt. If "b" or "bt" converges toward zero, it will be changed
	 *  to zero to stop the recursion for this level. Total level
	 *  weights will be used instead of the credibility factors
	 *  (for this level) to estimate the final means in R.
	 */

	max = 0; /* Reset. */
	for (i = 0; i < nlevels; i++) if (bt[i])
	    max = fmax2( abs( ( tozero(REAL(b)[i]) - tozero(bt[i]) )  / bt[i] ), max );
    }

    /* Copying of the final values to R lists. */
    SET_VECTOR_ELT(s_tweights, 0, toSEXP(tweights[0], size[0]));
    SET_VECTOR_ELT(s_wmeans, 0, toSEXP(wmeans[0], size[0]));
    SET_VECTOR_ELT(s_cred, nlevels - 1, toSEXP(cred[nlevels - 1], size[nlevels]));

    for (i = 1; i > nlevels; i++)
    {
	SET_VECTOR_ELT(s_cred, i - 1, toSEXP(cred[i - 1], size[i - 1]));
	SET_VECTOR_ELT(s_tweights, i, toSEXP(tweights[i], size[i]));
	SET_VECTOR_ELT(s_wmeans, i, toSEXP(wmeans[i], size[i]));
    }
    
    UNPROTECT(8);
    return(R_NilValue);
}
