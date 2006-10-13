/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Pareto distribution. See ../R/pareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dpareto(double x, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    return R_D_exp(log(shape) + shape * log(scale)
		   - (shape + 1.0) * log(x + scale));
}

double ppareto(double q, double shape, double scale, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    u = exp(-log1p(exp(log(limit) - log(scale))));

    return R_DT_Cval(R_pow(u, -shape));
}

double qpareto(double p, double shape, double scale, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale * (R_pow(R_D_Cval(p), -1.0 / shape) - 1.0);
}

double rpareto(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    return scale * (R_pow(unif_rand(), -1.0 / shape) - 1.0);
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
	* pbeta(u, tmp1, tmp2, 1, 0) / gammafn(shape)
	+ R_pow(limit, order) * R_pow(0.5 - u + 0.5, shape);
}
