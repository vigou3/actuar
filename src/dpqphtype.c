/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability density, cumulative probability
 *  and moment generating functions, and raw moments for phase-type
 *  distributions. This file is based on dpq.c with the following
 *  modifications:
 *
 *     1. support for a matrix argument;
 *     2. no iteration over the parameters;
 *     3. support for two parameter distributions only;
 *     4. many sanity checks on the arguments that are done in the
 *        {d,p,r,m,mgf} functions for other probability laws are done
 *        here because of item 2 above.
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

#define if_NA_dpqphtype2_set(y, x)                              \
    if      (ISNA (x) || naargs) y = NA_REAL;                   \
    else if (ISNAN(x) || nanargs) y = R_NaN;                    \
    else if (naflag) y = R_NaN;


static SEXP dpqphtype2_1(SEXP sx, SEXP sa, SEXP sb, SEXP sI, double (*f)())
{
    SEXP sy, bdims;
    int i, j, ij, n, m, sxo = OBJECT(sx);
    double tmp1, tmp2, *x, *a, *b, *y;
    int i_1;

    /* Flags used in sanity check of arguments. Listed from highest to
     * lowest priority. */
    Rboolean naargs = FALSE, nanargs = FALSE, naflag = FALSE;


#define SETUP_DPQPHTYPE2                                        \
    if (!isNumeric(sx) || !isNumeric(sa) || !isMatrix(sb))      \
        error(_("invalid arguments"));                          \
                                                                \
    n  = LENGTH(sx);                                            \
    if (n == 0)                                                 \
        return(allocVector(REALSXP, 0));                        \
                                                                \
    m = LENGTH(sa);                                             \
    bdims = getAttrib(sb, R_DimSymbol);                         \
    if (INTEGER(bdims)[0] != INTEGER(bdims)[1] ||               \
        INTEGER(bdims)[0] != m)                                 \
        naflag = TRUE;                                          \
                                                                \
    PROTECT(sx = coerceVector(sx, REALSXP));                    \
    PROTECT(sa = coerceVector(sa, REALSXP));                    \
    PROTECT(sb = coerceVector(sb, REALSXP));                    \
    PROTECT(sy = allocVector(REALSXP, n));                      \
    x = REAL(sx);                                               \
    a = REAL(sa);                                               \
    b = REAL(sb);                                               \
    y = REAL(sy);                                               \
                                                                \
    tmp1 = 0.0;                                                 \
    for (i = 0; i < m && !naargs && !nanargs && !naflag; i++)   \
    {                                                           \
        if ((naargs = ISNA(a[i])))                              \
            break;                                              \
        if ((nanargs = ISNAN(a[i])))                            \
            break;                                              \
        tmp1 += a[i];                                           \
        tmp2 = 0.0;                                             \
        for (j = 0; j < m; j++)                                 \
        {                                                       \
            ij = i + j * m;                                     \
            if ((naargs = ISNA(b[ij])))                         \
                break;                                          \
            if ((nanargs = ISNAN(b[ij])))                       \
                break;                                          \
            if (i == j && (naflag = b[ij] >= 0))                \
                break;                                          \
            if (i != j && (naflag = b[ij] < 0))                 \
                break;                                          \
            tmp2 += b[ij];                                      \
        }                                                       \
        if (!(naargs || nanargs))                               \
            naflag = tmp2 > 0;                                  \
    }                                                           \
    if (!(naargs || nanargs))                                   \
        naflag = tmp1 > 1


    SETUP_DPQPHTYPE2;

    i_1 = asInteger(sI);
    for (i = 0; i < n; i++)
    {
        if_NA_dpqphtype2_set(y[i], x[i])
        else
        {
            y[i] = f(x[i], a, b, m, i_1);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

#define FINISH_DPQPHTYPE2                               \
    if (naflag)                                         \
        warning(R_MSG_NA);                              \
                                                        \
    SET_ATTRIB(sy, duplicate(ATTRIB(sx)));              \
    SET_OBJECT(sy, sxo);                                \
                                                        \
    UNPROTECT(4)

    FINISH_DPQPHTYPE2;

    return sy;
}

static SEXP dpqphtype2_2(SEXP sx, SEXP sa, SEXP sb, SEXP sI, SEXP sJ, double (*f)())
{
    SEXP sy, bdims;
    int i, j, ij, n, m, sxo = OBJECT(sx);
    double tmp1, tmp2, *x, *a, *b, *y;
    int i_1, i_2;

    /* Flags used in sanity check of arguments. Listed from highest to
     * lowest priority. */
    Rboolean naargs = FALSE, nanargs = FALSE, naflag = FALSE;

    SETUP_DPQPHTYPE2;

    i_1 = asInteger(sI);
    i_2 = asInteger(sJ);
    for (i = 0; i < n; i++)
    {
        if_NA_dpqphtype2_set(y[i], x[i])
        else
        {
            y[i] = f(x[i], a, b, m, i_1, i_2);
            if (ISNAN(y[i])) naflag = TRUE;
        }
    }

    FINISH_DPQPHTYPE2;

    return sy;
}

#define DPQPHTYPE2_1(A, FUN) dpqphtype2_1(CAR(A), CADR(A), CADDR(A), CADDDR(A), FUN);
#define DPQPHTYPE2_2(A, FUN) dpqphtype2_2(CAR(A), CADR(A), CADDR(A), CADDDR(A), CAD4R(A), FUN)

SEXP actuar_do_dpqphtype2(int code, SEXP args)
{
    switch (code)
    {
    case  1:  return DPQPHTYPE2_1(args, dphtype);
    case  2:  return DPQPHTYPE2_2(args, pphtype);
    case  3:  return DPQPHTYPE2_1(args, mphtype);
    case  4:  return DPQPHTYPE2_1(args, mgfphtype);
    default:
        error(_("internal error in actuar_do_dpqphtype2"));
    }

    return args;                /* never used; to keep -Wall happy */
}

/* Main function, the only one used by .External(). */
SEXP actuar_do_dpqphtype(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to actuar_do_dpqphtype{1,2,3,4,5} */
    for (i = 0; dpq_tab[i].name; i++)
        if (!strcmp(dpq_tab[i].name, name))
            return dpq_tab[i].cfun(dpq_tab[i].code, CDR(args));

    /* No dispatch is an error */
    error("internal error in actuar_do_dpqphtype");

    return args;                /* never used; to keep -Wall happy */
}
