/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse Burr distribution. See ../R/invburr.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvburr(double x, double shape1, double shape2, double scale,
		int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape1 * shape2 * u^shape1 * (1 - u) / x
     *
     *  with u = v/(1 + v), v = (x/scale)^shape2.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp = shape2 * (log(x) - log(scale));
    log1mu = - log1p(exp(tmp));
    logu = tmp + log1mu;

    return R_D_exp(log(shape1) + log(shape2) + shape1 * logu
		   + log1mu - log(x));
}

double pinvburr(double q, double shape1, double shape2, double scale,
		int lower_tail, int log_p)
{
    double u, tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
	return 1;

    tmp = shape2 * (log(q) - log(scale));
    u = exp(tmp - log1p(exp(tmp)));

    return R_DT_val(R_pow(u, shape1));
}

double qinvburr(double p, double shape1, double shape2, double scale,
		int lower_tail, int log_p)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale * R_pow(R_pow(R_D_Lval(p), -1.0/shape1) - 1.0, -1.0/shape2);
}

double rinvburr(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    return scale * R_pow(R_pow(unif_rand(), -1.0/shape1) - 1.0, -1.0/shape2);
}

double minvburr(double order, double shape1, double shape2, double scale,
		int give_log)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= - shape1 * shape2 ||
	order  >= shape2)
	return R_NaN;

    tmp = order / shape2;

    return R_pow(scale, order) * gammafn(shape1 + tmp) * gammafn(1.0 - tmp)
	/ gammafn(shape1);
}

double levinvburr(double limit, double shape1, double shape2, double scale,
		  double order, int give_log)
{
    double u, tmp, tmp1, tmp2, tmp3;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= -shape1 * shape2)
	return R_NaN;

    if (limit <= 0.0)
	return 0;

    if (!R_FINITE(limit))
	return mburr(order, shape1, shape2, scale, 0);

    tmp1 = order / shape2;
    tmp2 = shape1 + tmp1;
    tmp3 = 1.0 - tmp1;

    tmp = shape2 * (log(limit) - log(scale));
    u = exp(tmp - log1p(exp(tmp)));

    return R_pow(scale, order) * gammafn(tmp2) * gammafn(tmp3)
	* pbeta(u, tmp2, tmp3, 1, 0) / gammafn(shape1)
	+ R_pow(limit, order) * (0.5 - R_pow(u, shape1) + 0.5);
}
