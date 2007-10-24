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

double dphtype(double x, double *alpha, double *S, int m, int give_log)
{
    /*  Density function is
     *
     *	alpha   * exp(x * S) * s
     *  (1 x m)   (m x m)      (m x 1)
     *
     *  with s = -S * e and e is a 1-vector.
     */

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    int i, j, jm;
    double *s, *Stmp;

    /* Build vector s (equal to minux the row sums of matrix S) and
     * matrix Stmp = x * S. Matrix S is stored in column major
     * order! */
    for (i = 0; i < m; i++)
    {
	s[i] = 0.0;
	for (j = 0; j < m; j++)
	{
	    jm = j * m;
	    s[i] -= S[i + jm];
	    Stmp[i + jm] = x * S[i + jm];
	}
    }

    return R_D_val(expmprod(alpha, Stmp, s, m));
}

double pphtype(double q, double *alpha, double *S, int m, int lower_tail,
	       int log_p)
{
    /*  Cumulative distribution function is
     *
     *	1 - alpha   * exp(q * S) * e
     *      (1 x m)   (m x m)      (m x 1)
     *
     *  where e is a 1-vector.
     */

    if (q <= 0)
	return R_DT_0;

    int i;
    double *e, *tmp;

    /* Create the 1-vector and multiply each element of S by q. */
    for (i = 0; i < m; i++)
	e[i] = 1;
    for (i = 0; i < m * m; i++)
	tmp[i] = q * S[i];

    return R_DT_Cval(expmprod(alpha, tmp, e, m));
}

double rphtype(double *alpha, double *S, int m)
{
    return 0.0;
}

double mphtype(double order, double *alpha, double *S, int m, int give_log)
{
    /*  Raw moment is
     *
     *	order!  * alpha   * (-S)^(-order) * 1
     *  (1 x 1)   (1 x m)   (m x m)         (m x 1)
     */

    if (order < 0.0 || 		/* negative moment */
	ftrunc(order) != order)	/* non-integer moment */
	return R_NaN;

    /* return R_D_val(matpow(S, m, (int) order)); */
    return order;
}

double mgfphtype(double x, double *alpha, double *S, int m, int give_log)
{
    return x;
}
