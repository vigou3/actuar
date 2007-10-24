/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute probability density, cumulative probability
 *  and moment generating functions, and raw moments for phase-type
 *  distributions. This file is based on dpq.c with the following
 *  modifications:
 *
 *     1. support for a matrix argument;
 *     2. no iteration over the parameters;
 *     3. support for two parameter distributions only.
 *
 *  Note that the "q" in the functions and file names was retained for
 *  symmetry reasons only, since the quantile function is not
 *  otherwise supported.
 *
 *  For details, see dpq.c.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"


static SEXP dpqphtype2_1(SEXP sx, SEXP sa, SEXP sb, SEXP sI, double (*f)())
{
    SEXP sy, bdims;
    int i, n, na, nb, nrow, ncol, sxo = OBJECT(sx);
    double *x, *a, *b, *y;
    int i_1;
    Rboolean naflag = FALSE;


#define SETUP_DPQPHTYPE2					\
    if (!isNumeric(sx) || !isNumeric(sa) || !isMatrix(sb))	\
	error(_("invalid arguments"));                          \
                                                                \
    bdims = getAttrib(sb, R_DimSymbol);                         \
    nrow = INTEGER(bdims)[0];                                   \
    ncol = INTEGER(bdims)[1];					\
                                                                \
    if (nrow != ncol)                                           \
	error(_("non-square transition matrix"));               \
								\
    n  = LENGTH(sx);						\
    na = LENGTH(sa);						\
    nb = LENGTH(sb);						\
								\
    if (na != nrow)                                             \
	error(_("non-conformable arguments"));			\
    if ((n == 0) || (na == 0) || (nb == 0))			\
	return(allocVector(REALSXP, 0));			\
								\
    PROTECT(sx = coerceVector(sx, REALSXP));			\
    PROTECT(sa = coerceVector(sa, REALSXP));			\
    PROTECT(sb = duplicate(sb));				\
    PROTECT(sy = allocVector(REALSXP, n));			\
    x = REAL(sx);						\
    a = REAL(sa);						\
    b = REAL(sb);						\
    y = REAL(sy);						\
								\
    if (na == 1)						\
	if (ISNA(a[0]) || ISNAN(a[0]))				\
	    naflag = TRUE;					\
    if (nb == 1)						\
	if (ISNA(b[0]) || ISNAN(b[0]))				\
	    naflag = TRUE

    SETUP_DPQPHTYPE2;

    if (naflag)
	for (i = 0; i < n; i++)
	    y[i] = a[0] + b[0];
    else
    {
	i_1 = asInteger(sI);
	for (i = 0; i < n; i++)
	{
	    y[i] = f(x[i], a, b, na, i_1);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

#define FINISH_DPQPHTYPE2				\
    if (naflag)						\
	warning(R_MSG_NA);				\
							\
    SET_ATTRIB(sy, duplicate(ATTRIB(sx)));		\
    SET_OBJECT(sy, sxo);				\
							\
    UNPROTECT(4)

    FINISH_DPQPHTYPE2;

    return sy;
}

static SEXP dpqphtype2_2(SEXP sx, SEXP sa, SEXP sb, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy, bdims;
    int i, n, na, nb, nrow, ncol, sxo = OBJECT(sx);
    double *x, *a, *b, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQPHTYPE2;

    if (naflag)
	for (i = 0; i < n; i++)
	    y[i] = a[0] + b[0];
    else
    {
	i_1 = asInteger(sI);
	i_2 = asInteger(sJ);
	for (i = 0; i < n; i++)
	{
	    y[i] = f(x[i], a, b, na, i_1, i_2);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

    FINISH_DPQPHTYPE2;

    return sy;
}

#define DPQPHTYPE2_1(A, FUN) dpqphtype2_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), FUN);
#define DPQPHTYPE2_2(A, FUN) dpqphtype2_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN)

SEXP do_dpqphtype2(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQPHTYPE2_1(args, dphtype);
    case  2:  return DPQPHTYPE2_2(args, pphtype);
    case  3:  return DPQPHTYPE2_1(args, mphtype);
    case  4:  return DPQPHTYPE2_1(args, mgfphtype);
    default:
	error(_("internal error in do_dpqphtype2"));
    }

    return args;		/* never used; to keep -Wall happy */
}

/* Main function, the only one used by .External(). */
SEXP do_dpqphtype(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to do_dpqphtype{1,2,3,4,5} */
    for (i = 0; fun_tab[i].name; i++)
    {
	if (!strcmp(fun_tab[i].name, name))
	{
	    return fun_tab[i].cfun(fun_tab[i].code, CDR(args));
	}
    }

    /* No dispatch is an error */
    error("internal error in do_dpqphtype");

    return args;		/* never used; to keep -Wall happy */
}
