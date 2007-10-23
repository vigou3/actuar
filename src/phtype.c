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
#include "locale.h"
#include "dpq.h"

double dphtype(double x, double *alpha, double *S, int m, int give_log)
{
    /*  Density function is
     *
     *	alpha   * exp(x * S) * s
     *  (1 x m)   (m x m)      (m x 1)
     *
     *  with s = S * e and e is a 1-vector.
     */

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    int i, j, im;
    double *s[m], *z;

    /* Build vector s --- equal to the row sums of matrix S --- and
     * multiply each element of S by x.  */
    for (i = 0; i < m; i++)
    {
	s[i] = 0.0;
	im = i * m;
	for (j = 0; j < m; j++)
	{
	    s[i] += S[im + j];
	    S[im + j] *= x;
	}
    }

    /* Compute alpha * exp(x * S) * s */
    expmprod(alpha, S, s, m, z)

    return R_D_val(z);
}

double pphtype(double x, double *alpha, double *S, int m, int lower_tail,
	       int log_p)
{
    /*  Cumulative distribution function is
     *
     *	1 - alpha   * exp(x * S) * e
     *      (1 x m)   (m x m)      (m x 1)
     *
     *  where e is a 1-vector.
     */

    if (q <= 0)
	return R_DT_0;

    int i;
    double *e[m], *z;

    /* Create the 1-vector and multiply each element of S by x. */
    for (i = 0; i < m; i++)
	e[i] = 1;
    for (i = 0; i < m * m; i++)
	S[i] *= x;

    /* Compute alpha * exp(x * S) * 1 */
    expmprod(alpha, S, e, m, z)

    return R_DT_Cval(z);
}

double rphtype(double *alpha, double *S, int m)
{
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

    int i, j, im;
    double *s[m], *z;

    /* Build vector s --- equal to the row sums of matrix S --- and
     * multiply each element of S by x.  */
    for (i = 0; i < m; i++)
    {
	s[i] = 0.0;
	im = i * m;
	for (j = 0; j < m; j++)
	{
	    s[i] += S[im + j];
	    S[im + j] *= x;
	}
    }

    /* Compute alpha * exp(x * S) * s */
    expmprod(alpha, S, s, m, z)

    return R_D_val(z);
}


double mpareto(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -1.0 ||
	order >= shape)
	return R_NaN;

    return R_pow(scale, order) * gammafn(1.0 + order) * gammafn(shape - order)
	/ gammafn(shape);
}

double levpareto(double limit, double shape, double scale, double order,
		 int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (limit <= 0.0)
	return 0;

    tmp1 = 1.0 + order;
    tmp2 = shape - order;

    u = exp(-log1p(exp(log(limit) - log(scale))));

    return R_pow(scale, order) * gammafn(tmp1) * gammafn(tmp2)
	* pbeta(0.5 - u + 0.5, tmp1, tmp2, 1, 0) / gammafn(shape)
	+ R_VG__0(limit, order) * R_pow(u, shape);
}
