/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to generate variates of some probability laws not in base
 *  R. Function .External() calls do_random() with arguments:
 *
 *       1. the name of the distribution from which to simulate, with
 *          an "r" prepended to it (e.g. "rpareto");
 *       2. the number of variates;
 *    3-... the parameters of the distribution.
 *
 *  Function do_random() will extract the name of the distribution,
 *  look up in table fun_tab defined in names.c which of
 *  do_random{1,2,3,4} should take care of the simulation and dispatch
 *  to this function. In turn, functions do_random{1,2,3,4} call
 *  function rdist() to get actual variates from distribution "dist".
 *
 *  This scheme is essentially what is used in base R (see files
 *  .../src/main/random.c, .../src/main/names.c in R sources). It looks
 *  convoluted at first, but adding a new distribution is plain simple:
 *  write an rdist() function, add an entry in names.c and in the
 *  definition of the corresponding do_random{1,2,3,4} function,
 *  declare the function in actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"

/* Functions for one parameter distributions */
static Rboolean random1(double (*f)(), double *a, int na, double *x, int n)
{
    double ai; 
    int i;
    Rboolean naflag = FALSE;
    for (i = 0; i < n; i++) 
    {
	ai = a[i % na];
	x[i] = f(ai);
	if (!R_FINITE(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RAND1(num, fun) \
	case num: \
	    random1(fun, REAL(a), na, REAL(x), n); \
	    break

SEXP do_random1(int code, SEXP args)
{
    SEXP x, a;
    int i, n, na;
    Rboolean naflag = FALSE;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) || !isNumeric(CADR(args)))
	error("invalid arguments");
    
    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1) 
    {
	n = asInteger(CAR(args));
	if (n == NA_INTEGER || n < 0)
	    error("invalid arguments");
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

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    if (na < 1)
    {
	for (i = 0; i < n; i++)
	    REAL(x)[i] = NA_REAL;
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	PROTECT(a = coerceVector(CADR(args), REALSXP));
	naflag = FALSE;
	GetRNGstate();
	switch (code) 
	{
	default:
	    error("internal error in do_random1");
	}
	if (naflag)
	    warning("NAs produced");
	PutRNGstate();
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for two parameter distributions */
static Rboolean random2(double (*f)(), double *a, int na, 
			double *b, int nb, double *x, int n)
{
    double ai, bi; 
    int i;
    Rboolean naflag = FALSE;
    for (i = 0; i < n; i++) 
    {
	ai = a[i % na];
	bi = b[i % nb];
	x[i] = f(ai, bi);
	if (!R_FINITE(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RAND2(num, fun) \
	case num: \
	    random2(fun, REAL(a), na, REAL(b), nb, REAL(x), n); \
	    break

SEXP do_random2(int code, SEXP args)
{
    SEXP x, a, b;
    int i, n, na, nb;
    Rboolean naflag = FALSE;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
	!isNumeric(CADR(args)) ||
	!isNumeric(CADDR(args)))
	error("invalid arguments");
    
    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1) 
    {
	n = asInteger(CAR(args));
	if (n == NA_INTEGER || n < 0)
	    error("invalid arguments");
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

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    if (na < 1 || nb < 1) 
    {
	for (i = 0; i < n; i++)
	    REAL(x)[i] = NA_REAL;
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	PROTECT(a = coerceVector(CADR(args), REALSXP));
	PROTECT(b = coerceVector(CADDR(args), REALSXP));
	naflag = FALSE;
	GetRNGstate();
	switch (code) 
	{
	    RAND2(0, rpareto);
	default:
	    error("internal error in do_random2");
	}
	if (naflag)
	    warning("NAs produced");
	PutRNGstate();
	UNPROTECT(2);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for three parameter distributions */
static Rboolean random3(double (*f) (), double *a, int na, 
			double *b, int nb, double *c, int nc, 
			double *x, int n)
{
    double ai, bi, ci;
    int i;
    Rboolean naflag = FALSE;
    for (i = 0; i < n; i++) 
    {
	ai = a[i % na];
	bi = b[i % nb];
	ci = c[i % nc];
	x[i] = f(ai, bi, ci);
	if (!R_FINITE(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RAND3(num, fun) \
	case num: \
	    random3(fun, REAL(a), na, REAL(b), nb, REAL(c), nc, REAL(x), n); \
	    break

SEXP do_random3(int code, SEXP args)
{
    SEXP x, a, b, c;
    int i, n, na, nb, nc;
    Rboolean naflag = FALSE;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
	!isNumeric(CADR(args)) ||
	!isNumeric(CADDR(args)) ||
	!isNumeric(CADDDR(args)))
	error("invalid arguments");
    
    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1) 
    {
	n = asInteger(CAR(args));
	if (n == NA_INTEGER || n < 0)
	    error("invalid arguments");
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

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    nc = LENGTH(CADDDR(args));
    if (na < 1 || nb < 1 || nc < 1) 
    {
	for (i = 0; i < n; i++)
	    REAL(x)[i] = NA_REAL;
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	PROTECT(a = coerceVector(CADR(args), REALSXP));
	PROTECT(b = coerceVector(CADDR(args), REALSXP));
	PROTECT(c = coerceVector(CADDDR(args), REALSXP));
	naflag = FALSE;
	GetRNGstate();
	switch (code) 
	{
	default:
	    error("internal error in do_random3");
	}
	if (naflag)
	    warning("NAs produced");
	PutRNGstate();
	UNPROTECT(3);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for four parameter distributions */
static Rboolean random4(double (*f) (), double *a, int na, 
			double *b, int nb, double *c, int nc, 
			double *d, int nd, double *x, int n)
{
    double ai, bi, ci, di;
    int i;
    Rboolean naflag = FALSE;
    for (i = 0; i < n; i++) 
    {
	ai = a[i % na];
	bi = b[i % nb];
	ci = c[i % nc];
	di = d[i % nd];
	x[i] = f(ai, bi, ci, di);
	if (!R_FINITE(x[i])) naflag = TRUE;
    }
    return(naflag);
}

#define RAND4(num, fun) \
	case num: \
	    random4(fun, REAL(a), na, REAL(b), nb, REAL(c), nc, REAL(d), nd, REAL(x), n); \
	    break

SEXP do_random4(int code, SEXP args)
{
    SEXP x, a, b, c, d;
    int i, n, na, nb, nc, nd;
    Rboolean naflag = FALSE;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
	!isNumeric(CADR(args)) ||
	!isNumeric(CADDR(args)) ||
	!isNumeric(CADDDR(args)) ||
	!isNumeric(CAD4R(args)))
	error("invalid arguments");
    
    /* Number of variates to generate */
    if (LENGTH(CAR(args)) == 1) 
    {
	n = asInteger(CAR(args));
	if (n == NA_INTEGER || n < 0)
	    error("invalid arguments");
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

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    nc = LENGTH(CADDDR(args));
    nd = LENGTH(CAD4R(args));
    if (na < 1 || nb < 1 || nc < 1 || nd < 1) 
    {
	for (i = 0; i < n; i++)
	    REAL(x)[i] = NA_REAL;
    }
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	PROTECT(a = coerceVector(CADR(args), REALSXP));
	PROTECT(b = coerceVector(CADDR(args), REALSXP));
	PROTECT(c = coerceVector(CADDDR(args), REALSXP));
	PROTECT(d = coerceVector(CAD4R(args), REALSXP));
	naflag = FALSE;
	GetRNGstate();
	switch (code) 
	{
	default:
	    error("internal error in do_random4");
	}
	if (naflag)
	    warning("NAs produced");
	PutRNGstate();
	UNPROTECT(4);
    }
    UNPROTECT(1);
    return x;
}



/* Main function, the only one used by .External(). */
SEXP do_random(SEXP args)
{
    int i;
    char *name;
    SEXP x;
    
    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to do_random{1,2,3,4} */
    for (i = 0; fun_tab[i].name; i++) 
    { 
	if (!strcmp(fun_tab[i].name, name))
	{
	    return fun_tab[i].cfun(fun_tab[i].code, CDR(args)); 
	}
    }

    /* If dispatch did not work, return NULL... */
    error("internal error in do_random");
    PROTECT(x = allocVector(NILSXP, 1));
    UNPROTECT(1);
    return(x);
}
