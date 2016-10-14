/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability density, cumulative probability
 *  quantile functions and moment generating functions, raw moments
 *  and limited moments for some probability laws not in base R (or
 *  those quantities not provided in base R). Function .External()
 *  calls actuar_do_dpq() with arguments:
 *
 *       1. the name of the distribution, with a "d", a "p" or "q"
 *          prepended to it (e.g. "dpareto", "pburr");
 *       2. the value(s) where the function is to be evaluated;
 *     3:x. the parameters of the distribution (including the order of
 *          the limited moment for lev*);
 *     x+1. whether to return the lower or upper tail probability or
 *          quantile (p* and q* only); see note below for m* and lev*
 *          functions;
 *     x+2. whether to return probability in log scale or the cumulant
 *          generating function (d*, p*, q* and mgf* only).
 *
 *  Function actuar_do_dpq() will extract the name of the distribution, look
 *  up in table dpq_tab defined in names.c which of actuar_do_dpq{1,2,3,4}
 *  should take care of the calculation and dispatch to this function.
 *  In turn, functions actuar_do_dpq{1,2,3,4} call function
 *  {d,p,q,m,lev,mgf}dist() to get actual values from distribution
 *  "dist".
 *
 *  Note: the m* and lev* functions came later in the process. In
 *  order to easily fit them into this system, I have decided to leave
 *  an unused 'give_log' argument in the C definitions of these
 *  functions. Otherwise, this would have required defining functions
 *  dpq{1,2,3,4,5}_0() below.
 *
 *  Functions therein are essentially identical to those found in
 *  .../src/main/arithmetic.c of R sources with a different naming
 *  scheme.
 *
 *  To add a new distribution: write a {d,p,q,m,lev,mgf}dist()
 *  function, add an entry in names.c and in the definition of the
 *  corresponding actuar_do_dpq{1,2,3,4} function, declare the function in
 *  actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much indirect help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"

/* Additional access macros */
#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))


/* Functions for special integrals with "zero parameter" */
#define if_NA_dpq0_set(y, x)		    \
        if      (ISNA (x)) y = NA_REAL;     \
        else if (ISNAN(x)) y = R_NaN;

static SEXP dpq0_1(SEXP sx, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, nx, sxo = OBJECT(sx);
    double xi, *x, *y;
    int i_1;
    Rboolean naflag = FALSE;

    if (!isNumeric(sx))
        error(_("invalid arguments"));

    nx = LENGTH(sx);
    if (nx == 0)
        return(allocVector(REALSXP, 0));
    PROTECT(sx = coerceVector(sx, REALSXP));
    PROTECT(sy = allocVector(REALSXP, nx));
    x = REAL(sx);
    y = REAL(sy);

    i_1 = asInteger(sI);

    for (i = 0; i < nx; i++)
    {
        xi = x[i];
        if_NA_dpq0_set(y[i], xi)
        else
        {
            y[i] = f(xi, i_1);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    if (naflag)
        warning(R_MSG_NA);

    SET_ATTRIB(sy, duplicate(ATTRIB(sx)));
    SET_OBJECT(sy, sxo);
    UNPROTECT(2);

    return sy;
}

#define DPQ0_1(A, FUN) dpq0_1(CAR(A), CADR(A), FUN);

SEXP actuar_do_dpq0(int code, SEXP args)
{
    switch (code)
    {
    case 201: return DPQ0_1(args, expint);   /* special integral */
    default:
        error(_("internal error in actuar_do_dpq0"));
    }

    return args;                /* never used; to keep -Wall happy */
}



/* Functions for one parameter distributions */
#define if_NA_dpq1_set(y, x, a)                         \
        if      (ISNA (x) || ISNA (a)) y = NA_REAL;     \
        else if (ISNAN(x) || ISNAN(a)) y = R_NaN;

#define mod_iterate1(n1, n2, i1, i2)            \
        for (i = i1 = i2 = 0; i < n;            \
             i1 = (++i1 == n1) ? 0 : i1,        \
             i2 = (++i2 == n2) ? 0 : i2,        \
             ++i)

static SEXP dpq1_1(SEXP sx, SEXP sa, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, n, nx, na, sxo = OBJECT(sx), sao = OBJECT(sa);
    double xi, ai, *x, *a, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ1                              \
    if (!isNumeric(sx) || !isNumeric(sa))       \
        error(_("invalid arguments"));          \
                                                \
    nx = LENGTH(sx);                            \
    na = LENGTH(sa);                            \
    if ((nx == 0) || (na == 0))                 \
        return(allocVector(REALSXP, 0));        \
    n = (nx < na) ? na : nx;                    \
    PROTECT(sx = coerceVector(sx, REALSXP));    \
    PROTECT(sa = coerceVector(sa, REALSXP));    \
    PROTECT(sy = allocVector(REALSXP, n));      \
    x = REAL(sx);                               \
    a = REAL(sa);                               \
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

#define FINISH_DPQ1                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
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

SEXP actuar_do_dpq1(int code, SEXP args)
{
    switch (code)
    {
    case   1: return DPQ1_1(args, mexp);
    case   2: return DPQ1_1(args, dinvexp);
    case   3: return DPQ1_2(args, pinvexp);
    case   4: return DPQ1_2(args, qinvexp);
    case   5: return DPQ1_1(args, minvexp);
    case   6: return DPQ1_1(args, mgfexp);
    case 101: return DPQ1_1(args, dlogarithmic);
    case 102: return DPQ1_2(args, plogarithmic);
    case 103: return DPQ1_2(args, qlogarithmic);
    case 104: return DPQ1_1(args, dztpois);
    case 105: return DPQ1_2(args, pztpois);
    case 106: return DPQ1_2(args, qztpois);
    case 107: return DPQ1_1(args, dztgeom);
    case 108: return DPQ1_2(args, pztgeom);
    case 109: return DPQ1_2(args, qztgeom);
    case 201: return DPQ1_1(args, gammaint); /* special integral */
    default:
        error(_("internal error in actuar_do_dpq1"));
    }

    return args;                /* never used; to keep -Wall happy */
}



/* Functions for two parameter distributions */
#define if_NA_dpq2_set(y, x, a, b)                              \
        if      (ISNA (x) || ISNA (a) || ISNA (b)) y = NA_REAL; \
        else if (ISNAN(x) || ISNAN(a) || ISNAN(b)) y = R_NaN;

#define mod_iterate2(n1, n2, n3, i1, i2, i3)    \
        for (i = i1 = i2 = i3 = 0; i < n;       \
             i1 = (++i1 == n1) ? 0 : i1,        \
             i2 = (++i2 == n2) ? 0 : i2,        \
             i3 = (++i3 == n3) ? 0 : i3,        \
             ++i)

static SEXP dpq2_1(SEXP sx, SEXP sa, SEXP sb, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, n, nx, na, nb,
        sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb);
    double xi, ai, bi, *x, *a, *b, *y;
    int i_1;
    Rboolean naflag = FALSE;


#define SETUP_DPQ2                                              \
    if (!isNumeric(sx) || !isNumeric(sa) || !isNumeric(sb))     \
        error(_("invalid arguments"));                          \
                                                                \
    nx = LENGTH(sx);                                            \
    na = LENGTH(sa);                                            \
    nb = LENGTH(sb);                                            \
    if ((nx == 0) || (na == 0) || (nb == 0))                    \
        return(allocVector(REALSXP, 0));                        \
    n = nx;                                                     \
    if (n < na) n = na;                                         \
    if (n < nb) n = nb;                                         \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
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

#define FINISH_DPQ2                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
    else if (n == nb) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));  \
        SET_OBJECT(sy, sbo);                    \
    }                                           \
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

/* This is needed for qinvgauss that has three additional parameters
 * for the tolerance, the maximum number of iterations and echoing of
 * the iterations. */
static SEXP dpq2_5(SEXP sx, SEXP sa, SEXP sb, SEXP sI, SEXP sJ,
		   SEXP sT, SEXP sM, SEXP sE, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, n, nx, na, nb,
        sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb);
    double xi, ai, bi, *x, *a, *b, *y;
    int i_1, i_2, i_4, i_5;
    double d_3;
    Rboolean naflag = FALSE;

    SETUP_DPQ2;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);
    d_3 = asReal(sT);
    i_4 = asInteger(sM);
    i_5 = asInteger(sE);

    mod_iterate2(nx, na, nb, ix, ia, ib)
    {
        xi = x[ix];
        ai = a[ia];
        bi = b[ib];
        if_NA_dpq2_set(y[i], xi, ai, bi)
        else
        {
            y[i] = f(xi, ai, bi, i_1, i_2, d_3, i_4, i_5);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    FINISH_DPQ2;

    return sy;
}

#define DPQ2_1(A, FUN) dpq2_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), FUN);
#define DPQ2_2(A, FUN) dpq2_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN)
#define DPQ2_5(A, FUN) dpq2_5(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), CAD6R(A), CAD7R(A), FUN)

SEXP actuar_do_dpq2(int code, SEXP args)
{
    switch (code)
    {
    case   1: return DPQ2_1(args, mgamma);
    case   2: return DPQ2_1(args, dinvgamma);
    case   3: return DPQ2_2(args, pinvgamma);
    case   4: return DPQ2_2(args, qinvgamma);
    case   5: return DPQ2_1(args, minvgamma);
    case   6: return DPQ2_1(args, dinvparalogis);
    case   7: return DPQ2_2(args, pinvparalogis);
    case   8: return DPQ2_2(args, qinvparalogis);
    case   9: return DPQ2_1(args, minvparalogis);
    case  10: return DPQ2_1(args, dinvpareto);
    case  11: return DPQ2_2(args, pinvpareto);
    case  12: return DPQ2_2(args, qinvpareto);
    case  13: return DPQ2_1(args, minvpareto);
    case  14: return DPQ2_1(args, dinvweibull);
    case  15: return DPQ2_2(args, pinvweibull);
    case  16: return DPQ2_2(args, qinvweibull);
    case  17: return DPQ2_1(args, minvweibull);
    case  18: return DPQ2_1(args, dlgamma);
    case  19: return DPQ2_2(args, plgamma);
    case  20: return DPQ2_2(args, qlgamma);
    case  21: return DPQ2_2(args, mlgamma);
    case  22: return DPQ2_1(args, dllogis);
    case  23: return DPQ2_2(args, pllogis);
    case  24: return DPQ2_2(args, qllogis);
    case  25: return DPQ2_1(args, mllogis);
    case  26: return DPQ2_1(args, mlnorm);
    case  27: return DPQ2_1(args, dparalogis);
    case  28: return DPQ2_2(args, pparalogis);
    case  29: return DPQ2_2(args, qparalogis);
    case  30: return DPQ2_1(args, mparalogis);
    case  31: return DPQ2_1(args, dpareto);
    case  32: return DPQ2_2(args, ppareto);
    case  33: return DPQ2_2(args, qpareto);
    case  34: return DPQ2_1(args, mpareto);
    case  35: return DPQ2_1(args, dpareto1);
    case  36: return DPQ2_2(args, ppareto1);
    case  37: return DPQ2_2(args, qpareto1);
    case  38: return DPQ2_1(args, mpareto1);
    case  39: return DPQ2_1(args, mweibull);
    case  40: return DPQ2_1(args, levexp);
    case  41: return DPQ2_1(args, levinvexp);
    case  42: return DPQ2_1(args, mbeta);
    case  43: return DPQ2_1(args, mgfgamma);
    case  44: return DPQ2_1(args, mgfnorm);
    case  45: return DPQ2_1(args, mgfunif);
    case  46: return DPQ2_1(args, mgfinvgamma);
    case  47: return DPQ2_1(args, mnorm);
    case  48: return DPQ2_1(args, mchisq);
    case  49: return DPQ2_1(args, mgfchisq);
    case  50: return DPQ2_1(args, minvGauss); /* deprecated v2.0-0 */
    case  51: return DPQ2_1(args, mgfinvGauss);
    case  52: return DPQ2_1(args, munif);
    case  53: return DPQ2_1(args, dgumbel);
    case  54: return DPQ2_2(args, pgumbel);
    case  55: return DPQ2_2(args, qgumbel);
    case  56: return DPQ2_1(args, mgumbel);
    case  57: return DPQ2_1(args, dinvgauss);
    case  58: return DPQ2_2(args, pinvgauss);
    case  59: return DPQ2_5(args, qinvgauss);
    case  60: return DPQ2_1(args, minvgauss);
    case 101: return DPQ2_1(args, dztnbinom);
    case 102: return DPQ2_2(args, pztnbinom);
    case 103: return DPQ2_2(args, qztnbinom);
    case 104: return DPQ2_1(args, dztbinom);
    case 105: return DPQ2_2(args, pztbinom);
    case 106: return DPQ2_2(args, qztbinom);
    case 107: return DPQ2_1(args, dzmlogarithmic);
    case 108: return DPQ2_2(args, pzmlogarithmic);
    case 109: return DPQ2_2(args, qzmlogarithmic);
    case 110: return DPQ2_1(args, dzmpois);
    case 111: return DPQ2_2(args, pzmpois);
    case 112: return DPQ2_2(args, qzmpois);
    case 113: return DPQ2_1(args, dzmgeom);
    case 114: return DPQ2_2(args, pzmgeom);
    case 115: return DPQ2_2(args, qzmgeom);
    case 201: return DPQ2_1(args, betaint); /* special integral */
    default:
        error(_("internal error in actuar_do_dpq2"));
    }

    return args;                /* never used; to keep -Wall happy */
}



/* Functions for three parameter distributions */
#define if_NA_dpq3_set(y, x, a, b, c)                                       \
        if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c)) y = NA_REAL; \
        else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c)) y = R_NaN;

#define mod_iterate3(n1, n2, n3, n4, i1, i2, i3, i4)    \
        for (i = i1 = i2 = i3 = i4 = 0; i < n;          \
             i1 = (++i1 == n1) ? 0 : i1,                \
             i2 = (++i2 == n2) ? 0 : i2,                \
             i3 = (++i3 == n3) ? 0 : i3,                \
             i4 = (++i4 == n4) ? 0 : i4,                \
             ++i)

static SEXP dpq3_1(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, n, nx, na, nb, nc,
        sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb), sco = OBJECT(sc);
    double xi, ai, bi, ci, *x, *a, *b, *c, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ3                                              \
    if (!isNumeric(sx) || !isNumeric(sa) ||                     \
        !isNumeric(sb) || !isNumeric(sc))                       \
        error(_("invalid arguments"));                          \
                                                                \
    nx = LENGTH(sx);                                            \
    na = LENGTH(sa);                                            \
    nb = LENGTH(sb);                                            \
    nc = LENGTH(sc);                                            \
    if ((nx == 0) || (na == 0) || (nb == 0) || (nc == 0))       \
        return(allocVector(REALSXP, 0));                        \
    n = nx;                                                     \
    if (n < na) n = na;                                         \
    if (n < nb) n = nb;                                         \
    if (n < nc) n = nc;                                         \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sc = coerceVector(sc, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
    c = REAL(sc);                                               \
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

#define FINISH_DPQ3                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
    else if (n == nb) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));  \
        SET_OBJECT(sy, sbo);                    \
    }                                           \
    else if (n == nc) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sc)));  \
        SET_OBJECT(sy, sco);                    \
    }                                           \
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

#define DPQ3_1(A, FUN) dpq3_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN);
#define DPQ3_2(A, FUN) dpq3_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), FUN)

SEXP actuar_do_dpq3(int code, SEXP args)
{
    switch (code)
    {
    case   1:  return DPQ3_1(args, dburr);
    case   2:  return DPQ3_2(args, pburr);
    case   3:  return DPQ3_2(args, qburr);
    case   4:  return DPQ3_1(args, mburr);
    case   5:  return DPQ3_1(args, dgenpareto);
    case   6:  return DPQ3_2(args, pgenpareto);
    case   7:  return DPQ3_2(args, qgenpareto);
    case   8:  return DPQ3_1(args, mgenpareto);
    case   9:  return DPQ3_1(args, dinvburr);
    case  10:  return DPQ3_2(args, pinvburr);
    case  11:  return DPQ3_2(args, qinvburr);
    case  12:  return DPQ3_1(args, minvburr);
    case  13:  return DPQ3_1(args, dinvtrgamma);
    case  14:  return DPQ3_2(args, pinvtrgamma);
    case  15:  return DPQ3_2(args, qinvtrgamma);
    case  16:  return DPQ3_1(args, minvtrgamma);
    case  17:  return DPQ3_1(args, dtrgamma);
    case  18:  return DPQ3_2(args, ptrgamma);
    case  19:  return DPQ3_2(args, qtrgamma);
    case  20:  return DPQ3_1(args, mtrgamma);
    case  21:  return DPQ3_1(args, levgamma);
    case  22:  return DPQ3_1(args, levinvgamma);
    case  23:  return DPQ3_1(args, levinvparalogis);
    case  24:  return DPQ3_1(args, levinvpareto);
    case  25:  return DPQ3_1(args, levinvweibull);
    case  26:  return DPQ3_1(args, levlgamma);
    case  27:  return DPQ3_1(args, levllogis);
    case  28:  return DPQ3_1(args, levlnorm);
    case  29:  return DPQ3_1(args, levparalogis);
    case  30:  return DPQ3_1(args, levpareto);
    case  31:  return DPQ3_1(args, levpareto1);
    case  32:  return DPQ3_1(args, levweibull);
    case  33:  return DPQ3_1(args, levbeta);
    case  34:  return DPQ3_1(args, levchisq);
    case  35:  return DPQ3_1(args, levinvGauss); /* notation from SuppDists */
    case  36:  return DPQ3_1(args, levunif);
    case 101:  return DPQ3_1(args, dzmnbinom);
    case 102:  return DPQ3_2(args, pzmnbinom);
    case 103:  return DPQ3_2(args, qzmnbinom);
    case 104:  return DPQ3_1(args, dzmbinom);
    case 105:  return DPQ3_2(args, pzmbinom);
    case 106:  return DPQ3_2(args, qzmbinom);
    default:
        error(_("internal error in actuar_do_dpq3"));
    }

    return args;                /* never used; to keep -Wall happy */
}



/* Functions for four parameter distributions */
#define if_NA_dpq4_set(y, x, a, b, c, d)                                   \
        if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c) || ISNA (d)) \
            y = NA_REAL;                                                   \
        else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c) || ISNAN(d)) \
            y = R_NaN;

#define mod_iterate4(n1, n2, n3, n4, n5, i1, i2, i3, i4, i5)    \
        for (i = i1 = i2 = i3 = i4 = i5 = 0; i < n;             \
             i1 = (++i1 == n1) ? 0 : i1,                        \
             i2 = (++i2 == n2) ? 0 : i2,                        \
             i3 = (++i3 == n3) ? 0 : i3,                        \
             i4 = (++i4 == n4) ? 0 : i4,                        \
             i5 = (++i5 == n5) ? 0 : i5,                        \
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

#define SETUP_DPQ4                                              \
    if (!isNumeric(sx) || !isNumeric(sa) || !isNumeric(sb) ||   \
        !isNumeric(sc) || !isNumeric(sd))                       \
        error(_("invalid arguments"));                          \
                                                                \
    nx = LENGTH(sx);                                            \
    na = LENGTH(sa);                                            \
    nb = LENGTH(sb);                                            \
    nc = LENGTH(sc);                                            \
    nd = LENGTH(sd);                                            \
    if ((nx == 0) || (na == 0) || (nb == 0) ||                  \
        (nc == 0) || (nd == 0))                                 \
        return(allocVector(REALSXP, 0));                        \
    n = nx;                                                     \
    if (n < na) n = na;                                         \
    if (n < nb) n = nb;                                         \
    if (n < nc) n = nc;                                         \
    if (n < nd) n = nd;                                         \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sc = coerceVector(sc, REALSXP));                    \
    PROTECT(sd = coerceVector(sd, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
    c = REAL(sc);                                               \
    d = REAL(sd);                                               \
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

#define FINISH_DPQ4                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
    else if (n == nb) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));  \
        SET_OBJECT(sy, sbo);                    \
    }                                           \
    else if (n == nc) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sc)));  \
        SET_OBJECT(sy, sco);                    \
    }                                           \
    else if (n == nd) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sd)));  \
        SET_OBJECT(sy, sdo);                    \
    }                                           \
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

#define DPQ4_1(A, FUN) dpq4_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), FUN);
#define DPQ4_2(A, FUN) dpq4_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), CAD6R(A), FUN)

SEXP actuar_do_dpq4(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ4_1(args, dtrbeta);
    case  2:  return DPQ4_2(args, ptrbeta);
    case  3:  return DPQ4_2(args, qtrbeta);
    case  4:  return DPQ4_1(args, mtrbeta);
    case  5:  return DPQ4_1(args, levburr);
    case  6:  return DPQ4_1(args, levgenpareto);
    case  7:  return DPQ4_1(args, levinvburr);
    case  8:  return DPQ4_1(args, levinvtrgamma);
    case  9:  return DPQ4_1(args, levtrgamma);
    case 10:  return DPQ4_1(args, dgenbeta);
    case 11:  return DPQ4_2(args, pgenbeta);
    case 12:  return DPQ4_2(args, qgenbeta);
    case 13:  return DPQ4_1(args, mgenbeta);
    default:
        error(_("internal error in actuar_do_dpq4"));
    }

    return args;                /* never used; to keep -Wall happy */
}

/* Functions for five parameter distributions */
#define if_NA_dpq5_set(y, x, a, b, c, d, e)                                \
        if      (ISNA (x) || ISNA (a) || ISNA (b) || ISNA (c) || ISNA (d) || ISNA (e)) \
            y = NA_REAL;                                                   \
        else if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c) || ISNAN(d) || ISNAN (e)) \
            y = R_NaN;

#define mod_iterate5(n1, n2, n3, n4, n5, n6, i1, i2, i3, i4, i5, i6)    \
        for (i = i1 = i2 = i3 = i4 = i5 = i6 = 0; i < n;                \
             i1 = (++i1 == n1) ? 0 : i1,                        \
             i2 = (++i2 == n2) ? 0 : i2,                        \
             i3 = (++i3 == n3) ? 0 : i3,                        \
             i4 = (++i4 == n4) ? 0 : i4,                        \
             i5 = (++i5 == n5) ? 0 : i5,                        \
             i6 = (++i6 == n6) ? 0 : i6,                        \
             ++i)

static SEXP dpq5_1(SEXP sx, SEXP sa, SEXP sb, SEXP sc, SEXP sd, SEXP se, SEXP sI, double (*f)())
{
    SEXP sy;
    int i, ix, ia, ib, ic, id, ie, n, nx, na, nb, nc, nd, ne,
        sxo = OBJECT(sx), sao = OBJECT(sa), sbo = OBJECT(sb),
        sco = OBJECT(sc), sdo = OBJECT(sd), seo = OBJECT(se);
    double xi, ai, bi, ci, di, ei, *x, *a, *b, *c, *d, *e, *y;
    int i_1;
    Rboolean naflag = FALSE;

#define SETUP_DPQ5                                              \
    if (!isNumeric(sx) || !isNumeric(sa) || !isNumeric(sb) ||   \
        !isNumeric(sc) || !isNumeric(sd) || !isNumeric(se))     \
        error(_("invalid arguments"));                          \
                                                                \
    nx = LENGTH(sx);                                            \
    na = LENGTH(sa);                                            \
    nb = LENGTH(sb);                                            \
    nc = LENGTH(sc);                                            \
    nd = LENGTH(sd);                                            \
    ne = LENGTH(se);                                            \
    if ((nx == 0) || (na == 0) || (nb == 0) ||                  \
        (nc == 0) || (nd == 0) || (ne == 0))                    \
        return(allocVector(REALSXP, 0));                        \
    n = nx;                                                     \
    if (n < na) n = na;                                         \
    if (n < nb) n = nb;                                         \
    if (n < nc) n = nc;                                         \
    if (n < nd) n = nd;                                         \
    if (n < ne) n = ne;                                         \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sc = coerceVector(sc, REALSXP));                    \
    PROTECT(sd = coerceVector(sd, REALSXP));                    \
    PROTECT(se = coerceVector(se, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
    c = REAL(sc);                                               \
    d = REAL(sd);                                               \
    e = REAL(se);                                               \
    y = REAL(sy)

    SETUP_DPQ5;

    i_1 = asInteger(sI);

    mod_iterate5(nx, na, nb, nc, nd, ne, ix, ia, ib, ic, id, ie)
    {
        xi = x[ix];
        ai = a[ia];
        bi = b[ib];
        ci = c[ic];
        di = d[id];
        ei = e[ie];
        if_NA_dpq5_set(y[i], xi, ai, bi, ci, di, ei)
        else
        {
            y[i] = f(xi, ai, bi, ci, di, ei, i_1);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

#define FINISH_DPQ5                             \
    if (naflag)                                 \
        warning(R_MSG_NA);                      \
                                                \
    if (n == nx) {                              \
        SET_ATTRIB(sy, duplicate(ATTRIB(sx)));  \
        SET_OBJECT(sy, sxo);                    \
    }                                           \
    else if (n == na) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sa)));  \
        SET_OBJECT(sy, sao);                    \
    }                                           \
    else if (n == nb) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sb)));  \
        SET_OBJECT(sy, sbo);                    \
    }                                           \
    else if (n == nc) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sc)));  \
        SET_OBJECT(sy, sco);                    \
    }                                           \
    else if (n == nd) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(sd)));  \
        SET_OBJECT(sy, sdo);                    \
    }                                           \
    else if (n == ne) {                         \
        SET_ATTRIB(sy, duplicate(ATTRIB(se)));  \
        SET_OBJECT(sy, seo);                    \
    }                                           \
    UNPROTECT(7)

    FINISH_DPQ5;

    return sy;
}

#define DPQ5_1(A, FUN) dpq5_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), CAD5R(A), CAD6R(A), FUN);

SEXP actuar_do_dpq5(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQ5_1(args, levtrbeta);
    case  2:  return DPQ5_1(args, levgenbeta);
    default:
        error(_("internal error in actuar_do_dpq5"));
    }

    return args;                /* never used; to keep -Wall happy */
}


/* Main function, the only one used by .External(). */
SEXP actuar_do_dpq(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to actuar_do_dpq{1,2,3,4,5} */
    for (i = 0; dpq_tab[i].name; i++)
    {
        if (!strcmp(dpq_tab[i].name, name))
        {
            return dpq_tab[i].cfun(dpq_tab[i].code, CDR(args));
        }
    }

    /* No dispatch is an error */
    error("internal error in actuar_do_dpq");

    return args;                /* never used; to keep -Wall happy */
}
