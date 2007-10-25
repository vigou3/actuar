/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and moment
 *  generating functions, raw moments and to simulate random variates
 *  for Phase-type distributions. See ../R/PhaseType.R for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "actuar.h"
#include "locale.h"
#include "dpq.h"

double dphtype(double x, double *pi, double *T, int m, int give_log)
{
    /*  Density function is
     *
     *	pi      * exp(x * T) * t
     *  (1 x m)   (m x m)      (m x 1)
     *
     *  with t = -T * e and e is a 1-vector.
     */

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    int i, j, jm;
    double *t, *tmp;

    /* Build vector t (equal to minux the row sums of matrix T) and
     * matrix tmp = x * T. */
    t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
    tmp = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	{
	    jm = j * m;
	    t[i] -= T[i + jm];
	    tmp[i + jm] = x * T[i + jm];
	}

    return R_D_val(expmprod(pi, tmp, t, m));
}

double pphtype(double q, double *pi, double *T, int m, int lower_tail,
	       int log_p)
{
    /*  Cumulative distribution function is
     *
     *	1 - pi      * exp(q * T) * e
     *      (1 x m)   (m x m)      (m x 1)
     *
     *  where e is a 1-vector.
     */

    if (q <= 0)
	return R_DT_0;

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

double rphtype(double *pi, double *T, int m)
{
    return 0.0;
}

double mphtype(double order, double *pi, double *T, int m, int give_log)
{
    /*  Raw moment is
     *
     *	order!  * pi      * (-T)^(-order) * e
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
    return x;
}
