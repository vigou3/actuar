/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute probability density, cumulative probability
 *  and quantile functions for some probability laws not in base
 *  R. Function .External() calls do_dpq() with arguments:
 *
 *       1. the name of the distribution from which to simulate, with
 *          a "d", a "p" or "q" prepended to it (e.g. "dpareto",
 *          "pburr");
 *       2. the value(s) where the function is to be evaluated;
 *     3-x. the parameters of the distribution;
 *     x+1. whether to return the lower or upper tail probability or
 *       quantile (p* and q* only);
 *     x+2. whether to return probability in log scale.
 *
 *  Function do_dpq() will extract the name of the distribution, look
 *  up in table fun_tab defined in names.c which of do_dpq{1,2,3,4}
 *  should take care of the calculation and dispatch to this
 *  function. In turn, functions do_dpq{1,2,3,4} call function
 *  {d,p,q}dist() to get actual values from distribution "dist".
 *
 *  Functions therein are essentially identical to those found in
 *  .../src/main/arithmetic.c of R sources with a different naming
 *  scheme.
 *
 *  To add a new distribution: write a {d,p,q}dist() function, add an
 *  entry in names.c and in the definition of the corresponding
 *  do_dpq{1,2,3,4} function, declare the function in actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"


/* Functions for two parameter distributions */
#define if_NA_dpq2_set(y, x, a, b)			        \
	if      (ISNA (x) || ISNA (a)|| ISNA (b)) y = NA_REAL;	\
	else if (ISNAN(x) || ISNAN(a)|| ISNAN(b)) y = R_NaN;

#define mod_iterate2(n1, n2, n3, i1, i2, i3)                    \
        for (i = i1 = i2 = i3 = 0; i < n;                       \
     	     i1 = (++i1 == n1) ? 0 : i1,			\
	     i2 = (++i2 == n2) ? 0 : i2,			\
	     i3 = (++i3 == n3) ? 0 : i3,			\
	     ++i)

static SEXP dpq2_1(SEXP sx, SEXP sa, SEXP sb, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, n, nx, na, nb,
	sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb);
    double xi, ai, bi, *x, *a, *b, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ2						\
    if (!isNumeric(sx) || !isNumeric(sa) || !isNumeric(sb))	\
	error(_("invalid arguments"));  			\
								\
    nx = LENGTH(sx);						\
    na = LENGTH(sa);						\
    nb = LENGTH(sb);						\
    if ((nx == 0) || (na == 0) || (nb == 0))			\
	return(allocVector(REALSXP, 0));			\
    n = nx;							\
    if (n < na) n = na;						\
    if (n < nb) n = nb;						\
    PROTECT(sx = coerceVector(sx, REALSXP));			\
    PROTECT(sa = coerceVector(sa, REALSXP));			\
    PROTECT(sb = coerceVector(sb, REALSXP));			\
    PROTECT(sy = allocVector(REALSXP, n));			\
    x = REAL(sx);						\
    a = REAL(sa);						\
    b = REAL(sb);						\
    y = REAL(sy)

    SETUP_DPQ2;

    i_1 = asInteger(sI);

    mod_iterate2(nx, na, nb, ix, ia, ib) {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
	if_NA_dpq2_set(y[i], xi, ai, bi)
	else 
	{
	    y[i] = f(xi, ai, bi, i_1);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

#define FINISH_DPQ2				\
    if(naflag)					\
	warning(_("NAs produced"));		\
						\
    if (n == nx) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sx)));	\
	SET_OBJECT(sy, sxo);			\
    }						\
    else if (n == na) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sa)));	\
	SET_OBJECT(sy, sao);			\
    }						\
    else if (n == nb) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sb)));	\
	SET_OBJECT(sy, sbo);			\
    }						\
    UNPROTECT(4)

    FINISH_DPQ2;

    return sy;
}

static SEXP dpq2_2(SEXP sx, SEXP sa, SEXP sb, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, n, nx, na, nb,
	sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb);
    double xi, ai, bi, *x, *a, *b, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQ2;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);

    for (i = 0; i < n; i++) 
    {
	xi = x[i % nx];
	ai = a[i % na];
	bi = b[i % nb];
	if_NA_dpq2_set(y[i], xi, ai, bi)
	else 
	{
	    y[i] = f(xi, ai, bi, i_1, i_2);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

    FINISH_DPQ2;

    return sy;
}

#define DPQ2_1(A, FUN) dpq2_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), FUN);
#define DPQ2_2(A, FUN) dpq2_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN)

SEXP do_dpq2(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ2_1(args, dpareto);
    case  2:  return DPQ2_2(args, ppareto);
    case  3:  return DPQ2_2(args, qpareto);
    case  4:  return DPQ2_1(args, dpareto1);
    case  5:  return DPQ2_2(args, ppareto1);
    case  6:  return DPQ2_2(args, qpareto1);
    default:
	error(_("internal error in do_dpq2"));
    }

    return args;		/* never used; to keep -Wall happy */
}


/* Main function, the only one used by .External(). */
SEXP do_dpq(SEXP args)
{
    int i;
    char *name;

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

    /* No dispatch is an error */
    error("internal error in do_random");

    return args;		/* never used; to keep -Wall happy */
}
