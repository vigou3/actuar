/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the paralogistic distribution. See ../R/paralogis.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dparalogis(double x, double shape, double scale, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp1 = shape * (log(x) - log(scale));
    tmp2 = R_pow(x / scale, shape);

    return  give_log ?
	2.0 * log(shape) + tmp1 - log(x) - (shape + 1.0) * log(1.0 + tmp2) :
	R_pow(shape, 2.0) * tmp2 / (x * R_pow(1.0 + tmp2, shape + 1.0));
}

double pparalogis(double q, double shape, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    tmp1 = R_pow(scale, shape);
    tmp2 = R_pow(q, shape);

    return lower_tail ?
	R_D_exp(log(1.0 - exp(shape * (shape * log(scale) - log(tmp1 + tmp2))))):
	R_D_exp(shape * (shape * log(scale) - log(tmp1 + tmp2)));
}

double qparalogis(double p, double shape, double scale, int lower_tail, int log_p)
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

    return lower_tail ? scale * R_pow(R_pow(1.0 - tmp, -tmp1) - 1.0, tmp1) :
	scale * R_pow(R_pow(tmp, -tmp1) - 1.0, tmp1);
}

double rparalogis(double shape, double scale)
{
    double a, tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    tmp = 1.0 / shape;

    a =  unif_rand();

    return scale * R_pow(R_pow(a, -tmp) - 1.0, tmp);
}

double mparalogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= R_pow(shape, 2.0))
	return R_NaN;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(shape - tmp) / gammafn(shape);
}

double levparalogis(double limit, double shape, double scale, double order, int give_log)
{
    double u, tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= R_pow(shape, 2.0) ||
	limit <= 0.0)
	return R_NaN;

    tmp = order / shape;
    u = 1.0 / (1.0 + R_pow(limit / scale, shape));

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(shape - tmp) * pbeta(1.0 - u, 1.0 + tmp, shape - tmp, 1, 0) + R_pow(limit, order) * R_pow(u, shape);
}
