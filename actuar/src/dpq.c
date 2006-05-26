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



/* Functions for one parameter distributions */
#define if_NA_dpq1_set(y, x, a)			        \
	if      (ISNA (x) || ISNA (a)) y = NA_REAL;	\
	else if (ISNAN(x) || ISNAN(a)) y = R_NaN;

#define mod_iterate1(n1, n2, i1, i2)		\
        for (i = i1 = i2 = 0; i < n;		\
     	     i1 = (++i1 == n1) ? 0 : i1,	\
	     i2 = (++i2 == n2) ? 0 : i2,	\
	     ++i)

static SEXP dpq1_1(SEXP sx, SEXP sa, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, n, nx, na, sxo = OBJECT(sx), sao = OBJECT(sa);
    double xi, ai, *x, *a, *b, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ1				\
    if (!isNumeric(sx) || !isNumeric(sa))	\
	error(_("invalid arguments"));  	\
						\
    nx = LENGTH(sx);				\
    na = LENGTH(sa);				\
    if ((nx == 0) || (na == 0))			\
	return(allocVector(REALSXP, 0));	\
    n = (nx < na) ? na : nx;			\
    PROTECT(sx = coerceVector(sx, REALSXP));	\
    PROTECT(sa = coerceVector(sa, REALSXP));	\
    PROTECT(sy = allocVector(REALSXP, n));	\
    x = REAL(sx);				\
    a = REAL(sa);				\
    y = REAL(sy)

    SETUP_DPQ1;

    i_1 = asInteger(sI);

    mod_iterate1(nx, na, ix, ia)
    {
	xi = x[ix];
	ai = a[ia];
	if_NA_dpq1_set(y[i], xi, ai)
	else
	{
	    y[i] = f(xi, ai, i_1);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

#define FINISH_DPQ1				\
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
    UNPROTECT(3)

    FINISH_DPQ1;

    return sy;
}

static SEXP dpq1_2(SEXP sx, SEXP sa, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy;
    int i, ix, ia, n, nx, na, sxo = OBJECT(sx), sao = OBJECT(sa);
    double xi, ai, *x, *a, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQ1;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);

    mod_iterate1(nx, na, ix, ia)
    {
	xi = x[ix];
	ai = a[ia];
	if_NA_dpq1_set(y[i], xi, ai)
	else
	{
	    y[i] = f(xi, ai, i_1, i_2);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

    FINISH_DPQ1;

    return sy;
}

#define DPQ1_1(A, FUN) dpq1_1(CAR(A), CADR(A), CADDR(A), FUN);
#define DPQ1_2(A, FUN) dpq1_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), FUN)

SEXP do_dpq1(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ1_1(args, diexp);
    case  2:  return DPQ1_2(args, piexp);
    case  3:  return DPQ1_2(args, qiexp);
    default:
	error(_("internal error in do_dpq1"));
    }

    return args;		/* never used; to keep -Wall happy */
}



/* Functions for two parameter distributions */
#define if_NA_dpq2_set(y, x, a, b)			        \
	if      (ISNA (x) || ISNA (a) || ISNA (b)) y = NA_REAL;	\
	else if (ISNAN(x) || ISNAN(a) || ISNAN(b)) y = R_NaN;

#define mod_iterate2(n1, n2, n3, i1, i2, i3)	\
        for (i = i1 = i2 = i3 = 0; i < n;	\
     	     i1 = (++i1 == n1) ? 0 : i1,	\
	     i2 = (++i2 == n2) ? 0 : i2,	\
	     i3 = (++i3 == n3) ? 0 : i3,	\
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

    mod_iterate2(nx, na, nb, ix, ia, ib)
    {
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

    mod_iterate2(nx, na, nb, ix, ia, ib)
    {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
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
    case  0:  return DPQ2_1(args, dpareto);
    case  1:  return DPQ2_2(args, ppareto);
    case  2:  return DPQ2_2(args, qpareto);
    case  3:  return DPQ2_1(args, dpareto1);
    case  4:  return DPQ2_2(args, ppareto1);
    case  5:  return DPQ2_2(args, qpareto1);
    default:
	error(_("internal error in do_dpq2"));
    }

    return args;		/* never used; to keep -Wall happy */
}



/* Functions for three parameter distributions */
#define if_NA_dpq3_set(y, x, a, b, c)			 	      	    \
	if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c)) y = NA_REAL; \
	else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c)) y = R_NaN;

#define mod_iterate3(n1, n2, n3, n4, i1, i2, i3, i4)	\
        for (i = i1 = i2 = i3 = i4 = 0; i < n;		\
     	     i1 = (++i1 == n1) ? 0 : i1,		\
	     i2 = (++i2 == n2) ? 0 : i2,		\
	     i3 = (++i3 == n3) ? 0 : i3,		\
	     i4 = (++i4 == n4) ? 0 : i4,		\
	     ++i)

static SEXP dpq3_1(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, n, nx, na, nb, nc,
	sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb), sco = OBJECT(sc);
    double xi, ai, bi, ci, *x, *a, *b, *c, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ3						\
    if (!isNumeric(sx) || !isNumeric(sa) || 			\
	!isNumeric(sb) || !isNumeric(sc))			\
	error(_("invalid arguments"));  			\
								\
    nx = LENGTH(sx);						\
    na = LENGTH(sa);						\
    nb = LENGTH(sb);						\
    nc = LENGTH(sc);						\
    if ((nx == 0) || (na == 0) || (nb == 0) || (nc == 0))	\
	return(allocVector(REALSXP, 0));			\
    n = nx;							\
    if (n < na) n = na;						\
    if (n < nb) n = nb;						\
    if (n < nc) n = nc;						\
    PROTECT(sx = coerceVector(sx, REALSXP));			\
    PROTECT(sa = coerceVector(sa, REALSXP));			\
    PROTECT(sb = coerceVector(sb, REALSXP));			\
    PROTECT(sc = coerceVector(sc, REALSXP));			\
    PROTECT(sy = allocVector(REALSXP, n));			\
    x = REAL(sx);						\
    a = REAL(sa);						\
    b = REAL(sb);						\
    c = REAL(sc);						\
    y = REAL(sy)

    SETUP_DPQ3;

    i_1 = asInteger(sI);

    mod_iterate3(nx, na, nb, nc, ix, ia, ib, ic)
    {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
	ci = c[ic];
	if_NA_dpq3_set(y[i], xi, ai, bi, ci)
	else
	{
	    y[i] = f(xi, ai, bi, ci, i_1);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

#define FINISH_DPQ3				\
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
    else if (n == nc) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sc)));	\
	SET_OBJECT(sy, sco);			\
    }						\
    UNPROTECT(5)

    FINISH_DPQ3;

    return sy;
}

static SEXP dpq3_2(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, n, nx, na, nb, nc,
	sxo = OBJECT(sx), sao = OBJECT(sa),
	sbo = OBJECT(sb), sco = OBJECT(sc);
    double xi, ai, bi, ci, *x, *a, *b, *c, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQ3;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);

    mod_iterate3(nx, na, nb, nc, ix, ia, ib, ic)
    {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
	ci = c[ic];
	if_NA_dpq3_set(y[i], xi, ai, bi, ci)
	else
	{
	    y[i] = f(xi, ai, bi, ci, i_1, i_2);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

    FINISH_DPQ3;

    return sy;
}

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define DPQ3_1(A, FUN) dpq3_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN);
#define DPQ3_2(A, FUN) dpq3_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), FUN)

SEXP do_dpq3(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ3_1(args, dburr);
    case  2:  return DPQ3_2(args, pburr);
    case  3:  return DPQ3_2(args, qburr);
    default:
	error(_("internal error in do_dpq3"));
    }

    return args;		/* never used; to keep -Wall happy */
}



/* Functions for four parameter distributions */
#define if_NA_dpq4_set(y, x, a, b, c, d)				   \
	if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c) || ISNA (d)) \
            y = NA_REAL; 						   \
	else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c) || ISNAN(d)) \
	    y = R_NaN;

#define mod_iterate4(n1, n2, n3, n4, n5, i1, i2, i3, i4, i5)	\
        for (i = i1 = i2 = i3 = i4 = i5 = 0; i < n;		\
     	     i1 = (++i1 == n1) ? 0 : i1,			\
	     i2 = (++i2 == n2) ? 0 : i2,			\
	     i3 = (++i3 == n3) ? 0 : i3,			\
	     i4 = (++i4 == n4) ? 0 : i4,			\
	     i5 = (++i5 == n5) ? 0 : i5,			\
	     ++i)

static SEXP dpq4_1(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sd, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, id, n, nx, na, nb, nc, nd,
	sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb),
	sco = OBJECT(sc), sdo = OBJECT(sd);
    double xi, ai, bi, ci, di, *x, *a, *b, *c, *d, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ4						\
    if (!isNumeric(sx) || !isNumeric(sa) || !isNumeric(sb) ||	\
	!isNumeric(sc) || !isNumeric(sd))			\
	error(_("invalid arguments"));  			\
								\
    nx = LENGTH(sx);						\
    na = LENGTH(sa);						\
    nb = LENGTH(sb);						\
    nc = LENGTH(sc);						\
    nd = LENGTH(sd);						\
    if ((nx == 0) || (na == 0) || (nb == 0) || 			\
	(nc == 0) || (nd == 0))					\
	return(allocVector(REALSXP, 0));			\
    n = nx;							\
    if (n < na) n = na;						\
    if (n < nb) n = nb;						\
    if (n < nc) n = nc;						\
    if (n < nd) n = nd;						\
    PROTECT(sx = coerceVector(sx, REALSXP));			\
    PROTECT(sa = coerceVector(sa, REALSXP));			\
    PROTECT(sb = coerceVector(sb, REALSXP));			\
    PROTECT(sc = coerceVector(sc, REALSXP));			\
    PROTECT(sd = coerceVector(sd, REALSXP));			\
    PROTECT(sy = allocVector(REALSXP, n));			\
    x = REAL(sx);						\
    a = REAL(sa);						\
    b = REAL(sb);						\
    c = REAL(sc);						\
    d = REAL(sd);						\
    y = REAL(sy)

    SETUP_DPQ4;

    i_1 = asInteger(sI);

    mod_iterate4(nx, na, nb, nc, nd, ix, ia, ib, ic, id)
    {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
	ci = c[ic];
	di = d[id];
	if_NA_dpq4_set(y[i], xi, ai, bi, ci, di)
	else
	{
	    y[i] = f(xi, ai, bi, ci, di, i_1);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

#define FINISH_DPQ4				\
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
    else if (n == nc) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sc)));	\
	SET_OBJECT(sy, sco);			\
    }						\
    else if (n == nd) {				\
	SET_ATTRIB(sy, duplicate(ATTRIB(sd)));	\
	SET_OBJECT(sy, sdo);			\
    }						\
    UNPROTECT(6)

    FINISH_DPQ4;

    return sy;
}

static SEXP dpq4_2(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sd, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, id, n, nx, na, nb, nc, nd,
	sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb),
	sco = OBJECT(sc), sdo = OBJECT(sd);
    double xi, ai, bi, ci, di, *x, *a, *b, *c, *d, *y;
    int i_1, i_2;
    Rboolean naflag = FALSE;

    SETUP_DPQ4;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);

    mod_iterate4(nx, na, nb, nc, nd, ix, ia, ib, ic, id)
    {
	xi = x[ix];
	ai = a[ia];
	bi = b[ib];
	ci = c[ic];
	di = d[id];
	if_NA_dpq4_set(y[i], xi, ai, bi, ci, di)
	else
	{
	    y[i] = f(xi, ai, bi, ci, di, i_1, i_2);
	    if (ISNAN(y[i])) naflag = TRUE;
	}
    }

    FINISH_DPQ4;

    return sy;
}

#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define DPQ4_1(A, FUN) dpq4_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), FUN);
#define DPQ4_2(A, FUN) dpq4_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), CAD6R(A), FUN)

SEXP do_dpq4(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ4_1(args, dtrbeta);
    case  2:  return DPQ4_2(args, ptrbeta);
    case  3:  return DPQ4_2(args, qtrbeta);
    default:
	error(_("internal error in do_dpq3"));
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
