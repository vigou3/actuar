/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to generate variates of some probability laws not in
 *  base R. Function .External() calls actuar_do_random() with
 *  arguments:
 *
 *      1. the name of the distribution from which to simulate, with
 *         an "r" prepended to it (e.g. "rpareto");
 *      2. the number of variates;
 *    3:x. the parameters of the distribution.
 *
 *  Function actuar_do_random() will extract the name of the
 *  distribution, look up in table random_tab defined in names.c which of
 *  actuar_do_random{1,2,3,4} should take care of the simulation and
 *  dispatch to this function. In turn, functions
 *  actuar_do_random{1,2,3,4} call function rdist() to get actual
 *  variates from distribution "dist".
 *
 *  This scheme is essentially what is used in base R (see files
 *  .../src/main/random.c, .../src/main/names.c in R sources) with
 *  add-ons taken from .../src/library/stats/src/random.c to support
 *  return values that can be either real or integer.
 *
 * To add a new distribution: write an rdist() function, add an entry
 * in names.c and in the definition of the corresponding
 * actuar_do_random{1,2,3,4} function, declare the function in
 * actuar.h.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much indirect help from the R Core Team
 */

#include <R.h>
#include <Rinternals.h>
#include "actuar.h"
#include "locale.h"

/* Utility function used in actuar_do_random{1,2,3,4}. */
static void fill_with_NAs(SEXP x, int n, SEXPTYPE type) {
    int i;

    if (type == INTSXP) {
        for (i = 0; i < n; i++) {
            INTEGER(x)[i] = NA_INTEGER;
        }
    } else { /* REALSXP */
        for (i = 0; i < n; i++) {
            REAL(x)[i] = NA_REAL;
        }
    }
    warning(_("NAs produced"));
}


/* Functions for one parameter distributions */
static Rboolean random1(double (*f)(), double *a, int na, SEXP x, int n, SEXPTYPE type)
{
    int i;
    Rboolean naflag = FALSE;
    if (type == INTSXP)
    {
	double rx;
	int *ix = INTEGER(x);

	for (i = 0; i < n; i++)
	{
	    rx = f(a[i % na]);
	    if (ISNAN(rx) || rx > INT_MAX || rx <= INT_MIN)
	    {
		naflag = TRUE;
		ix[i] = NA_INTEGER;
	    }
	    else
		ix[i] = (int) rx;
	}
    }
    else /* REALSXP */
    {
	double *rx = REAL(x);
	for (i = 0; i < n; i++)
	{
	    rx[i] = f(a[i % na]);
	    if (ISNAN(rx[i])) naflag = TRUE;
	}
    }

    return(naflag);
}

#define RAND1(num, fun)                                 \
    case num:						\
    naflag = random1(fun, REAL(a), na, x, n, type);	\
    break

SEXP actuar_do_random1(int code, SEXP args, SEXPTYPE type)
{
    SEXP x, a;
    int n, na;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) || !isNumeric(CADR(args)))
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
    PROTECT(x = allocVector(type, n));
    if (n == 0)
    {
        UNPROTECT(1);
        return(x);
    }

    /* If length of parameters < 1, return NaN */
    na = LENGTH(CADR(args));
    if (na < 1)
	fill_with_NAs(x, n, type);
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
    	Rboolean naflag = FALSE;
        PROTECT(a = coerceVector(CADR(args), REALSXP));
        GetRNGstate();
        switch (code)
        {
            RAND1(1, rinvexp);
            RAND1(2, rlogarithmic);
            RAND1(3, rztpois);
            RAND1(4, rztgeom);
        default:
            error(_("internal error in actuar_do_random1"));
        }
        if (naflag)
            warning(R_MSG_NA);
        PutRNGstate();
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for two parameter distributions */
static Rboolean random2(double (*f)(), double *a, int na,
                        double *b, int nb, SEXP x, int n, SEXPTYPE type)
{
    int i;
    Rboolean naflag = FALSE;
    if (type == INTSXP)
    {
	double rx;
	int *ix = INTEGER(x);

	for (i = 0; i < n; i++)
	{
	    rx = f(a[i % na], b[i % nb]);
	    if (ISNAN(rx) || rx > INT_MAX || rx <= INT_MIN)
	    {
		naflag = TRUE;
		ix[i] = NA_INTEGER;
	    }
	    else
		ix[i] = (int) rx;
	}
    }
    else /* REALSXP */
    {
	double *rx = REAL(x);
	for (i = 0; i < n; i++)
	{
	    rx[i] = f(a[i % na], b[i % nb]);
	    if (ISNAN(rx[i])) naflag = TRUE;
	}
    }

    return(naflag);
}

#define RAND2(num, fun)							\
    case num:								\
    naflag = random2(fun, REAL(a), na, REAL(b), nb, x, n, type);	\
    break

SEXP actuar_do_random2(int code, SEXP args, SEXPTYPE type)
{
    SEXP x, a, b;
    int n, na, nb;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
        !isNumeric(CADR(args)) ||
        !isNumeric(CADDR(args)))
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
    PROTECT(x = allocVector(type, n));
    if (n == 0)
    {
        UNPROTECT(1);
        return(x);
    }

    /* If length of parameters < 1, return NA */
    na = LENGTH(CADR(args));
    nb = LENGTH(CADDR(args));
    if (na < 1 || nb < 1)
	fill_with_NAs(x, n, type);
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
    	Rboolean naflag = FALSE;
        PROTECT(a = coerceVector(CADR(args), REALSXP));
        PROTECT(b = coerceVector(CADDR(args), REALSXP));
        GetRNGstate();
        switch (code)
        {
            RAND2( 1, rinvgamma);
            RAND2( 2, rinvparalogis);
            RAND2( 3, rinvpareto);
            RAND2( 4, rinvweibull);
            RAND2( 5, rlgamma);
            RAND2( 6, rllogis);
            RAND2( 7, rparalogis);
            RAND2( 8, rpareto);
            RAND2( 9, rpareto1);
            RAND2(10, rztnbinom);
            RAND2(11, rztbinom);
            RAND2(12, rzmlogarithmic);
            RAND2(13, rzmpois);
        default:
            error(_("internal error in actuar_do_random2"));
        }
        if (naflag)
            warning(R_MSG_NA);
        PutRNGstate();
        UNPROTECT(2);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for three parameter distributions */
static Rboolean random3(double (*f) (), double *a, int na,
                        double *b, int nb, double *c, int nc,
                        SEXP x, int n, SEXPTYPE type)
{
    int i;
    Rboolean naflag = FALSE;
    if (type == INTSXP)
    {
	double rx;
	int *ix = INTEGER(x);

	for (i = 0; i < n; i++)
	{
	    rx = f(a[i % na], b[i % nb], c[i % nc]);
	    if (ISNAN(rx) || rx > INT_MAX || rx <= INT_MIN)
	    {
		naflag = TRUE;
		ix[i] = NA_INTEGER;
	    }
	    else
		ix[i] = (int) rx;
	}
    }
    else /* REALSXP */
    {
	double *rx = REAL(x);
	for (i = 0; i < n; i++)
	{
	    rx[i] = f(a[i % na], b[i % nb], c[i % nc]);
	    if (ISNAN(rx[i])) naflag = TRUE;
	}
    }

    return(naflag);
}

#define RAND3(num, fun)							\
    case num:								\
    naflag = random3(fun, REAL(a), na, REAL(b), nb, REAL(c), nc, x, n, type); \
    break

SEXP actuar_do_random3(int code, SEXP args, SEXPTYPE type)
{
    SEXP x, a, b, c;
    int n, na, nb, nc;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
        !isNumeric(CADR(args)) ||
        !isNumeric(CADDR(args)) ||
        !isNumeric(CADDDR(args)))
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
    PROTECT(x = allocVector(type, n));
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
	fill_with_NAs(x, n, type);
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
    	Rboolean naflag = FALSE;
        PROTECT(a = coerceVector(CADR(args), REALSXP));
        PROTECT(b = coerceVector(CADDR(args), REALSXP));
        PROTECT(c = coerceVector(CADDDR(args), REALSXP));
        GetRNGstate();
        switch (code)
        {
            RAND3(1, rburr);
            RAND3(2, rgenpareto);
            RAND3(3, rinvburr);
            RAND3(4, rinvtrgamma);
            RAND3(5, rtrgamma);
            RAND3(6, rzmnbinom);
            RAND3(7, rzmbinom);
        default:
            error(_("internal error in actuar_do_random3"));
        }
        if (naflag)
            warning(R_MSG_NA);
        PutRNGstate();
        UNPROTECT(3);
    }
    UNPROTECT(1);
    return x;
}


/* Functions for four parameter distributions */
static Rboolean random4(double (*f) (), double *a, int na,
                        double *b, int nb, double *c, int nc,
                        double *d, int nd, SEXP x, int n, SEXPTYPE type)
{
    int i;
    Rboolean naflag = FALSE;
    if (type == INTSXP)
    {
	double rx;
	int *ix = INTEGER(x);

	for (i = 0; i < n; i++)
	{
	    rx = f(a[i % na], b[i % nb], c[i % nc], d[i % nd]);
	    if (ISNAN(rx) || rx > INT_MAX || rx <= INT_MIN)
	    {
		naflag = TRUE;
		ix[i] = NA_INTEGER;
	    }
	    else
		ix[i] = (int) rx;
	}
    }
    else /* REALSXP */
    {
	double *rx = REAL(x);
	for (i = 0; i < n; i++)
	{
	    rx[i] = f(a[i % na], b[i % nb], c[i % nc], d[i % nd]);
	    if (ISNAN(rx[i])) naflag = TRUE;
	}
    }

    return(naflag);
}

#define RAND4(num, fun)							\
    case num:								\
    naflag = random4(fun, REAL(a), na, REAL(b), nb, REAL(c), nc, REAL(d), nd, x, n, type); \
    break

SEXP actuar_do_random4(int code, SEXP args, SEXPTYPE type)
{
    SEXP x, a, b, c, d;
    int n, na, nb, nc, nd;

    /* Check validity of arguments */
    if (!isVector(CAR(args)) ||
        !isNumeric(CADR(args)) ||
        !isNumeric(CADDR(args)) ||
        !isNumeric(CADDDR(args)) ||
        !isNumeric(CAD4R(args)))
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
    PROTECT(x = allocVector(type, n));
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
	fill_with_NAs(x, n, type);
    /* Otherwise, dispatch to appropriate r* function */
    else
    {
	Rboolean naflag = FALSE;
        PROTECT(a = coerceVector(CADR(args), REALSXP));
        PROTECT(b = coerceVector(CADDR(args), REALSXP));
        PROTECT(c = coerceVector(CADDDR(args), REALSXP));
        PROTECT(d = coerceVector(CAD4R(args), REALSXP));
        GetRNGstate();
        switch (code)
        {
            RAND4(1, rtrbeta);
            RAND4(2, rgenbeta);
        default:
            error(_("internal error in actuar_do_random4"));
        }
        if (naflag)
            warning(R_MSG_NA);
        PutRNGstate();
        UNPROTECT(4);
    }
    UNPROTECT(1);
    return x;
}



/* Main function, the only one used by .External(). */
SEXP actuar_do_random(SEXP args)
{
    int i;
    const char *name;

    /* Extract distribution name */
    args = CDR(args);
    name = CHAR(STRING_ELT(CAR(args), 0));

    /* Dispatch to actuar_do_random{1,2,3,4} */
    for (i = 0; random_tab[i].name; i++)
    {
        if (!strcmp(random_tab[i].name, name))
            return random_tab[i].cfun(random_tab[i].code, CDR(args), random_tab[i].type);
    }

    /* No dispatch is an error */
    error(_("internal error in actuar_do_random"));

    return args;                /* never used; to keep -Wall happy */
}
