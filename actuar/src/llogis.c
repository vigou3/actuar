/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the loglogistic distribution. See ../R/llogis.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dllogis(double x, double shape, double scale, int give_log)
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

    return give_log ?
	log(shape) + tmp1 - log(x) - 2.0 * log(1.0 + exp(tmp1)) :
	shape * tmp2 / (x * R_pow(1.0 + tmp2, 2.0));
}

double pllogis(double q, double shape, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
	return 1;

    tmp1 = R_pow(q / scale, shape);
    tmp2 = shape * (log(q) - log(scale));

    return lower_tail ? R_D_exp(tmp1 - log(1.0 + tmp2)):
	R_D_exp(log(1.0 - tmp1 / (1.0 + tmp1)));
}

double qllogis(double p, double shape, double scale, int lower_tail, int log_p)
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

    return lower_tail ? scale * R_pow(tmp / (1.0 - tmp), tmp1) :
	scale * R_pow((1.0 - tmp) / tmp, tmp1);
}

double rllogis(double shape, double scale)
{
    double a;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    a = unif_rand();

    return scale * R_pow(a / (1.0 - a), 1.0 / shape);
}

double mllogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= shape)
	return R_NaN;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(1.0 - tmp);
}

double levllogis(double limit, double shape, double scale, double order, int give_log)
{
    double tmp1, tmp2, tmp3;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(x) ||
	!R_FINITE(limit) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= shape ||
	limit <= 0.0)
	return R_NaN;

    tmp1 = R_pow(limit / scale, shape);
    tmp2 = tmp1 / (1.0 + tmp1);
    tmp3 = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp3) * gammafn(1.0 - tmp3) * pbeta(tmp2, 1.0 + tmp3, 1.0 - tmp3, 1, 0) + R_pow(limit, order) * (1.0 - tmp2);
}
