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

    return give_log ?
	log(shape) + shape * log(scale) - (shape + 1.0) * log(x + scale) :
	shape * R_pow(scale, shape) / R_pow(x + scale, shape + 1.0);
}

double ppareto(double x, double shape, double scale, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (x <= 0)
	return R_DT_0;

    tmp = log(scale) - log(x + scale);

    return lower_tail ? R_D_exp(log(1.0 - exp(shape * (tmp)))) :
	R_D_exp(shape * (tmp));
}

double qpareto(double p, double shape, double scale, int lower_tail, int log_p)
{
    double tmp, tmp1;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);
    tmp1 = 1.0 / shape;

    return lower_tail ? scale * (R_pow(1.0 - tmp, -tmp1) - 1.0) :
	scale * (R_pow(tmp, -tmp1) - 1.0);
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

    return R_pow(scale, order) * gammafn(order + 1.0) * gammafn(shape - order) / gammafn(shape);
}

double levpareto(double limit, double shape, double scale, double order, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	limit <= 0.0 ||
	order <= -1.0 ||
	order >= shape)
	return R_NaN;

    return R_pow(scale, order) * gammafn(order + 1.0) * gammafn(shape - order) * pbeta(limit / (limit + scale), order + 1.0, shape - order, 1, 0) / gammafn(shape) + R_pow(limit, order) * R_pow(scale / (scale + limit), shape);
}
