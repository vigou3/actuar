/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to generate variates of phase-type distributions. This
 *  file is based on random.c with the following modifications:
 *
 *     1. support for a matrix argument;
 *     2. no iteration over the parameters;
 *     3. support for two parameter distributions only.
 *
 *  For details, see random.c.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include "actuar.h"
#include "locale.h"

static Rboolean randomphtype2(double (*f)(), double *a, double *b,
			      int na, double *x, int n)
{
    int i, j;
    double *rates, **Q;
    Rboolean naflag = FALSE;

    /* The sub-intensity matrix and initial probability vector never
     * change, so compute the transition matrix of the Markov chain
     * and the vector of rate parameters before looping. */
    rates = (double *) R_alloc(na, sizeof(double));
    Q = (double *) R_alloc(na, sizeof(double));
    for (i = 0; i < na; i++)
    {
	Q[i] = (double *) S_alloc(na, sizeof(double));
	rates[i] = -b[i * (na + 1)];
	for (j = 0; j < na; j++)
	    if (i != j)
		Q[i][j] = b[i + j * na] / rates[i];
    }

    for (i = 0; i < n; i++)
    {
	x[i] = f(a, Q, rates, na);
	if (!R_FINITE(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RANDPHTYPE2(num, fun) \
	case num: \
	    randomphtype2(fun, REAL(a), REAL(b), na, REAL(x), n); \
	    break

SEXP do_randomphtype2(int code, SEXP args)
{
    SEXP x, a, b, bdims;
    int i, n, na, nb, nrow, ncol;
    Rboolean naflag = FALSE;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
	!isNumeric(CADR(args)) ||
	!isMatrix(CADDR(args)))
	error(_("invalid arguments"));

    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1)
    {
	n = asInteger(CAR(args));
	if (n == NA_INTEGER || n < 0)
	    error(_("invalid arguments"));
    }
    else
	n = LENGTH(CAR(args));

    /* If n == 0, return numeric(0) */
    PROTECT(x = allocVector(REALSXP, n));
    if (n == 0)
    {
	UNPROTECT(1);
	return(x);
    }

    /* Sanity checks of arguments. */
    PROTECT(a = coerceVector(CADR(args), REALSXP));
    PROTECT(b = coerceVector(CADDR(args), REALSXP));
    bdims = getAttrib(b, R_DimSymbol);
    nrow = INTEGER(bdims)[0];
    ncol = INTEGER(bdims)[1];
    if (nrow != ncol)
	error(_("non-square transition matrix"));
    na = LENGTH(a);
    nb = LENGTH(b);
    if (na != nrow)
	error(_("non-conformable arguments"));

    /*  If length of parameters < 1, or either of the two parameters
     *  is NA return NA. */
    if (na < 1 ||
	(na == 1 && !(R_FINITE(REAL(a)[0]) && R_FINITE(REAL(b)[0]))))
    {
	for (i = 0; i < n; i++)
	    REAL(x)[i] = NA_REAL;
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	naflag = FALSE;
	GetRNGstate();
	switch (code)
	{
	    RANDPHTYPE2(1, rphtype);
	default:
	    error(_("internal error in do_randomphtype2"));
	}
	if (naflag)
	    warning(R_MSG_NA);
	PutRNGstate();
    }
    UNPROTECT(3);
    return x;
}


/* Main function, the only one used by .External(). */
SEXP do_randomphtype(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to do_random{1,2,3,4} */
    for (i = 0; fun_tab[i].name; i++)
	if (!strcmp(fun_tab[i].name, name))
	    return fun_tab[i].cfun(fun_tab[i].code, CDR(args));

    /* No dispatch is an error */
    error(_("internal error in do_randomphtype"));

    return args;		/* never used; to keep -Wall happy */
}
