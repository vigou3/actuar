/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Function to compute the iterative part of function cm, used
 *  to deal with credibility models.
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

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define weights(i, j) (cred[i][j] != 0 ? cred[i][j] : tweights[i + 1][j])

SEXP toSEXP(double *x, int size)
{
    SEXP ans = allocVector(REALSXP, size);
    memcpy(REAL(ans), x, size * sizeof(double));
    return ans;
}

SEXP actuar_do_hierarc(SEXP args)
{
    SEXP s_cred, s_tweights, s_wmeans, s_fnodes, denoms, b, tol, maxit, echo;
    double **cred, **tweights, **wmeans, diff, bw;
    int **fnodes, nlevels, i, j, k, count = 0;

    /*  All values received from R are protected. */
    PROTECT(s_cred = coerceVector(CADR(args), VECSXP));
    PROTECT(s_tweights = coerceVector(CADDR(args), VECSXP));
    PROTECT(s_wmeans = coerceVector(CADDDR(args), VECSXP));
    PROTECT(s_fnodes = coerceVector(CAD4R(args), VECSXP));
    PROTECT(denoms = coerceVector(CAD5R(args), REALSXP));
    PROTECT(b = coerceVector(CAD6R(args), REALSXP));
    PROTECT(tol = coerceVector(CAD7R(args), REALSXP));
    PROTECT(maxit = coerceVector(CAD8R(args), INTSXP));
    PROTECT(echo = coerceVector(CAD9R(args), LGLSXP));

    /* Initialization of some variables */
    double bt[length(b)];       /* previous values of 'b' */
    nlevels = length(b) - 1;    /* number of levels in the model */
    bt[nlevels] = REAL(b)[nlevels]; /* within entity variance; never
                                     * changes. */
    int size[nlevels + 1];      /* total number of nodes at each
                                 * level, including the portfolio level */
    size[0] = 1;
    for (i = 1; i <= nlevels; i++)
        size[i] = length(VECTOR_ELT(s_fnodes, i - 1));

    /* Allocation of arrays that will be needed below. */
    cred     = (double **) R_alloc(nlevels,     sizeof(double *));
    tweights = (double **) R_alloc(nlevels + 1, sizeof(double *));
    wmeans   = (double **) R_alloc(nlevels + 1, sizeof(double *));
    fnodes   = (int **)    R_alloc(nlevels,     sizeof(int *));
    tweights[0] = (double *) R_alloc(size[0], sizeof(double));
    wmeans[0]   = (double *) R_alloc(size[0], sizeof(double));
    for (i = 1; i <= nlevels; i++)
    {
        cred[i - 1]   = (double *) R_alloc(size[i], sizeof(double));
        tweights[i]   = (double *) R_alloc(size[i], sizeof(double));
        wmeans[i]     = (double *) R_alloc(size[i], sizeof(double));
        fnodes[i - 1] = (int *)    R_alloc(size[i], sizeof(int));
    }

    /* Get values of fnodes, tweights and wmeans from R lists. For
     * the latter two, only the entity level values are initialized
     * in R or meaningful. */
    for (i = 0; i < nlevels; i++)
        memcpy(fnodes[i], INTEGER(VECTOR_ELT(s_fnodes, i)),
               size[i + 1] * sizeof(int));
    memcpy(tweights[nlevels], REAL(VECTOR_ELT(s_tweights, nlevels)),
           size[nlevels] * sizeof(double));
    memcpy(wmeans[nlevels], REAL(VECTOR_ELT(s_wmeans, nlevels)),
           size[nlevels] * sizeof(double));

    /* If printing of iterations was asked for, start by printing a
     * header and the starting values. */
    if (LOGICAL(echo)[0])
    {
        Rprintf("Iteration\tVariance estimates\n %d\t\t", count);
        for (i = 0; i < nlevels; i++)
            Rprintf(" %.8g  ", REAL(b)[i]);
        Rprintf("\n");
    }

    /* Iterative part. */
    do
    {
        /* Stop after 'maxit' iterations and issue warning. */
        if (++count > INTEGER(maxit)[0])
        {
            warning(_("maximum number of iterations reached before obtaining convergence"));
            break;
        }

        /* Copy the previous values of 'b'. */
        for (i = 0; i < nlevels; i++)
            bt[i] = REAL(b)[i];

        /* Run through all levels from lowest to highest. */
        for (i = nlevels - 1; i >= 0; i--)
        {
            /* Reset the total weights and weighted averages. */
            for (j = 0; j < size[i]; j++)
            {
                tweights[i][j] = 0;
                wmeans[i][j] = 0;
            }

            /* Find the first non-zero within variance estimator. */
            for (j = 1; REAL(b)[i + j] == 0; j++);
            bw = REAL(b)[i + j];

            /* Calculation of the new credibility factors, total
             * weights and (numerators of) weighted averages. */
            for (j = 0; j < size[i + 1]; j++)
            {
                cred[i][j] = 1 / (1 + bw / (REAL(b)[i] * tweights[i + 1][j]));
                k = fnodes[i][j] - 1; /* C version of tapply(). */
                tweights[i][k] += weights(i, j);
                wmeans[i][k] += weights(i, j) * wmeans[i + 1][j];
            }

            /* Final calculation of weighted averages with the
             * division by the total weight. */
            for (j = 0; j < size[i]; j++)
            {
                if (tweights[i][j] > 0)
                    wmeans[i][j] = wmeans[i][j] / tweights[i][j];
                else
                    wmeans[i][j] = 0;
            }

            /* Calculation of the new current level variance estimator
             * only if the previous one is strictly positive. */
            if (bt[i] > 0)
            {
                REAL(b)[i] = 0;
                for (j = 0; j < size[i + 1]; j++)
                {
                    k = fnodes[i][j];
                    REAL(b)[i] += weights(i, j) * R_pow_di(wmeans[i + 1][j] - wmeans[i][k - 1], 2);
                }
                REAL(b)[i] = REAL(b)[i] / REAL(denoms)[i];

                /* Set the estimator to 0 if it is close enough to 0
                 * and henceforth stop iterations on this
                 * parameter. */
                if (REAL(b)[i] <= R_pow_di(REAL(tol)[0], 2))
                    REAL(b)[i] = 0;
            }

            /* Recompute the credibility factors, total weights and
             * weighted means with the latest between variance
             * estimator. */
            for (j = 0; j < size[i]; j++)
            {
                tweights[i][j] = 0;
                wmeans[i][j] = 0;
            }
            for (j = 0; j < size[i + 1]; j++)
            {
                cred[i][j] = 1 / (1 + bw / (REAL(b)[i] * tweights[i + 1][j]));
                k = fnodes[i][j] - 1;
                tweights[i][k] += weights(i, j);
                wmeans[i][k] += weights(i, j) * wmeans[i + 1][j];
            }
            for (j = 0; j < size[i]; j++)
            {
                if (tweights[i][j] > 0)
                    wmeans[i][j] = wmeans[i][j] / tweights[i][j];
                else
                    wmeans[i][j] = 0;
            }
        }

        /* Trace */
        if (LOGICAL(echo)[0])
        {
            Rprintf(" %d\t\t", count);
            for (i = 0; i < nlevels; i++)
                Rprintf(" %.8g  ", REAL(b)[i]);
            Rprintf("\n");
        }

        /*  Computation of the largest difference between two
         *  iterations. Estimators set to 0 are not taken into
         *  account. */
        diff = 0;
        for (i = 0; i < nlevels; i++)
            if (REAL(b)[i] > 0)
                diff = fmax2(abs(REAL(b)[i] - bt[i])/bt[i], diff);
    }
    while (diff >= REAL(tol)[0]);

    /* Copy the final values to R lists. */
    SET_VECTOR_ELT(s_tweights, 0, toSEXP(tweights[0], size[0]));
    SET_VECTOR_ELT(s_wmeans,   0, toSEXP(wmeans[0], size[0]));
    for (i = 1; i <= nlevels; i++)
    {
        SET_VECTOR_ELT(s_cred,     i - 1, toSEXP(cred[i - 1], size[i]));
        SET_VECTOR_ELT(s_tweights, i,     toSEXP(tweights[i], size[i]));
        SET_VECTOR_ELT(s_wmeans,   i,     toSEXP(wmeans[i],   size[i]));
    }

    UNPROTECT(9);
    return(R_NilValue);
}
