/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and moment
 *  generating functions, raw moments and to simulate random variates
 *  for Phase-type distributions. See ../R/PhaseType.R for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include "actuar.h"
#include "locale.h"
#include "dpq.h"

double dphtype(double x, double *pi, double *T, int m, int give_log)
{
    /*  Density function is
     *
     *  pi      * exp(x * T) * t
     *  (1 x m)   (m x m)      (m x 1)
     *
     *  for x > 0, with t = -T * e and e a 1-vector, and 1 - pi * e
     *  for x = 0.
     */

    if (!R_FINITE(x) || x < 0.0)
        return R_D__0;

    if (x == 0.0)
    {
        int i;
        double z = 0.0;

        for (i = 0; i < m; i++)
            z += pi[i];

        return R_D_Clog(z);
    }

    int i, j, ij;
    double *t, *tmp;

    /* Build vector t (equal to minus the row sums of matrix T) and
     * matrix tmp = x * T. */
    t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
    tmp = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
        {
            ij = i + j * m;
            t[i] -= T[ij];
            tmp[ij] = x * T[ij];
        }

    return R_D_val(expmprod(pi, tmp, t, m));
}

double pphtype(double q, double *pi, double *T, int m, int lower_tail,
               int log_p)
{
    /*  Cumulative distribution function is
     *
     *  1 - pi      * exp(q * T) * e
     *      (1 x m)   (m x m)      (m x 1)
     *
     *  for x > 0, where e a 1-vector, and 1 - pi * e for x = 0.
     */

    if (q < 0.0)
        return R_DT_0;

    if (q == 0.0)
    {
        int i;
        double z = 0.0;

        for (i = 0; i < m; i++)
            z += pi[i];

        return R_DT_Cval(z);
    }

    int i;
    double *e, *tmp;

    /* Create the 1-vector and multiply each element of T by q. */
    e = (double *) R_alloc(m, sizeof(double));
    for (i = 0; i < m; i++)
        e[i] = 1;
    tmp = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m * m; i++)
        tmp[i] = q * T[i];

    return R_DT_Cval(expmprod(pi, tmp, e, m));
}

double rphtype(double *pi, double **Q, double *rates, int m)
{
    /* Algorithm based on Neuts, M. F. (1981), "Generating random
     * variates from a distribution of phase type", WSC '81:
     * Proceedings of the 13th conference on Winter simulation, IEEE
     * Press, <http://portal.acm.org/citation.cfm?id=802607&coll=portal&dl=ACM#>
     */

    int i, j, state, *nvisits;
    double z = 0.0;

    nvisits = (int *) S_alloc(m, sizeof(int));

    /* Simulate initial state according to vector pi (transient states
     * are numbered 0, ..., m - 1 and absorbing state is numbered
     * m). See the definition of SampleSingleValue() to see why this
     * works fine here and below. */
    state = SampleSingleValue(m, pi);

    /* Simulate the underlying Markov chain using transition matrix Q
     * while counting the number of visits in each transient state. */
    while (state != m)
    {
        nvisits[state]++;
        state = SampleSingleValue(m, Q[state]);
    }

    /* Variate is the sum of as many exponential variates as there are
     * visits in each state, with the rate parameter varying per
     * state. */
    for (i = 0; i < m; i++)
        for (j = 0; j < nvisits[i]; j++)
            z += exp_rand() / rates[i];

    return z;
}

double mphtype(double order, double *pi, double *T, int m, int give_log)
{
    /*  Raw moment is
     *
     *  order!  * pi      * (-T)^(-order) * e
     *  (1 x 1)   (1 x m)   (m x m)         (m x 1)
     *
     * where e is a 1-vector. Below, the moment is computed as
     * (-1)^order * order! * sum(pi * T^(-order))
     */

    if (order < 0.0 || (int) order != order)
        return R_NaN;

    int i, j;
    double tmp = 0.0, *Tpow;

    /* Compute the power of T */
    Tpow = (double *) R_alloc(m * m, sizeof(double));
    matpow(T, m, (int) -order, Tpow);

    /* Compute vector tmp = sum(pi * Tpow) */
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
            tmp += pi[j] * Tpow[i * m + j];

    /* Multiply by -1 if order is odd */
    return R_D_val((int) order % 2 ?
                   -gammafn(order + 1.0) * tmp :
                   gammafn(order + 1.0) * tmp);
}

double mgfphtype(double x, double *pi, double *T, int m, int give_log)
{
    /*  Moment generating function is
     *
     *  pi      * (-x * I - T)^(-1) * t       + (1 - pi      * e)
     *  (1 x m)   (m x m)             (m x 1)        (1 x m)   (m x 1)
     *
     *  with t = -T * e, e a 1-vector and I the identity matrix.
     *  Below, the mgf is computed as 1 - pi * (e + (x * I + T)^(-1) * t.
     */

    if (x == 0.0)
        return R_D_exp(0.0);

    int i, j, ij;
    double z = 0.0, *t, *tmp1, *tmp2;

    /* Build vector t (equal to minux the row sums of matrix T) and
     * matrix tmp1 = x * I + T. */
    t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
    tmp1 = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m; i++)
        for (j = 0; j < m; j++)
        {
            ij = i + j * m;
            t[i] -= T[ij];
            tmp1[ij] = (i == j) ? x + T[ij] : T[ij];
        }

    /* Compute tmp2 = tmp1^(-1) * t */
    tmp2 = (double *) R_alloc(m, sizeof(double));
    solve(tmp1, t, m, 1, tmp2);

    /* Compute z = pi * (e + tmp2) */
    for (i = 0; i < m; i++)
        z += pi[i] * (1 + tmp2[i]);

    return R_D_Clog(z);
}
