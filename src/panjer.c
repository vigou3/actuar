/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Function to compute the recursive part of the Panjer formula
 *  to approximate the aggregate claim amount distribution of
 *  a portfolio over a period.
 *
 *  AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "locale.h"

#define CAD5R(e) CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))
#define CAD8R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))
#define CAD9R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))))

SEXP do_panjer(SEXP args)
{
    SEXP p0, p1, fs0, sfx, a, b, tol, maxit, echo, sfs;
    double *fs, *fx, cumul;
    int upper, m, k, x = 1;
    double norm;                /* normalizing constant */
    double term;                /* constant in the (a, b, 1) case */

    /*  The length of vector fs is not known in advance. We opt for a
     *  simple scheme: allocate memory for a vector of size 'size',
     *  double the size when the vector is full. */
    int size = 100;
    fs = (double *) S_alloc(size, sizeof(double));
    for (k = 0; k < size; k++)
        fs[k] = 0;

    /*  All values received from R are then protected. */
    PROTECT(p0 = coerceVector(CADR(args), REALSXP));
    PROTECT(p1 = coerceVector(CADDR(args), REALSXP));
    PROTECT(fs0 = coerceVector(CADDDR(args), REALSXP));
    PROTECT(sfx = coerceVector(CAD4R(args), REALSXP));
    PROTECT(a = coerceVector(CAD5R(args), REALSXP));
    PROTECT(b = coerceVector(CAD6R(args), REALSXP));
    PROTECT(tol = coerceVector(CAD7R(args), REALSXP));
    PROTECT(maxit = coerceVector(CAD8R(args), INTSXP));
    PROTECT(echo = coerceVector(CAD9R(args), LGLSXP));

    /* Initialization of some variables */
    fx = REAL(sfx);             /* severity distribution */
    upper = length(sfx) - 1;    /* upper bound of the support of the
                                 * severity distribution */
    fs[0] = REAL(fs0)[0];       /* value of Pr[S = 0] (computed in R) */
    cumul = REAL(fs0)[0];       /* cumulative probability computed */
    norm = 1 - REAL(a)[0] * fx[0]; /* normalizing constant */

    /* If printing of recursions was asked for, start by printing a
     * header and the probability at 0. */
    if (LOGICAL(echo)[0])
        Rprintf("x\tPr[S = x]\tCumulative probability\n%d\t%.8g\t%.8g\n",
                0, fs[0], fs[0]);

    /* (a, b, 0) case (if p0 is NULL) */
    if (isNull(CADR(args)))
        do
        {
            /* Stop after 'maxit' recursions and issue warning. */
            if (x > INTEGER(maxit)[0])
            {
                warning(_("maximum number of recursions reached before the probability distribution was complete"));
                break;
            }

            /* If fs is too small, double its size */
            if (x >= size)
            {
                fs = (double *) S_realloc((char *) fs, 2 * size, size, sizeof(double));
                for (k = size; k < 2 * size; k++)
                    fs[k] = 0;
                size = 2 * size;
            }

            m = fmin2(x, upper); /* upper bound of the sum */

            /* Compute probability up to the scaling constant */
            for (k = 1; k <= m; k++)
                fs[x] += (REAL(a)[0] + REAL(b)[0] * k / x) * fx[k] * fs[x - k];
            fs[x] = fs[x] / norm; /* normalization */
            cumul += fs[x];       /* cumulative sum */

            if (LOGICAL(echo)[0])
                Rprintf("%d\t%.8g\t%.8g\n", x, fs[x], cumul);

            x++;
        } while (cumul < REAL(tol)[0]);
    /* (a, b, 1) case (if p0 is non-NULL) */
    else
    {
        /* Next line is a hack to reproduce the fact that the
         * distribution of claim amounts is 0 past its maximum value.
         * Only needed in the (a, b, 1) case for the additional term
         * in the recursion formula. */
        fx[++upper] = 0;

        /* Constant term in the (a, b, 1) case. */
        term = (REAL(p1)[0] - (REAL(a)[0] + REAL(b)[0]) * REAL(p0)[0]);

        do
        {
            /* Stop after 'maxit' recursions and issue warning. */
            if (x > INTEGER(maxit)[0])
            {
                warning(_("maximum number of recursions reached before the probability distribution was complete"));
                break;
            }

            if (x >= size)
            {
                fs = (double *) S_realloc((char *) fs, 2 * size, size, sizeof(double));
                for (k = size; k < 2 * size; k++)
                    fs[k] = 0;
                size = 2 * size;
            }

            m = fmin2(x, upper);

            for (k = 1; k <= m; k++)
                fs[x] += (REAL(a)[0] + REAL(b)[0] * k / x) * fx[k] * fs[x - k];
            fs[x] = (fs[x] + fx[m] * term) / norm;
            cumul += fs[x];

            if (LOGICAL(echo)[0])
                Rprintf("%d\t%.8g\t%.8g\n", x, fs[x], cumul);

            x++;
        }
        while (cumul < REAL(tol)[0]);
    }

    /*  Copy the values of fs to a SEXP which will be returned to R. */
    PROTECT(sfs = allocVector(REALSXP, x));
    memcpy(REAL(sfs), fs, x * sizeof(double));

    UNPROTECT(10);
    return(sfs);
}
