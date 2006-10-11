/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse Weibull distribution. See ../R/invweibull.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvweibull(double x, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(shape)  ||
	scale <= 0.0 ||
	shape <= 0.0)
	return R_NaN;;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp = R_pow(scale / x, shape);

    return  give_log ?
	log(shape) + shape * (log(scale) - log(x)) - tmp - log(x) :
	shape * tmp * exp(-tmp) / x;
}

double pinvweibull(double q, double shape, double scale, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    tmp = R_pow(scale / q, shape);

    return lower_tail ? R_D_exp(-tmp): R_D_exp(log(1.0 - exp(-tmp)));
}

double qinvweibull(double p, double shape, double scale, int lower_tail, int log_p)
{
    double tmp, tmp1;

    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);

    tmp1 = -1.0 / shape;

    return lower_tail ?
	scale * R_pow(-log(tmp), tmp1):
	scale * R_pow(-log(1.0 - tmp), tmp1);
}

double rinvweibull(double shape, double scale)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	return R_NaN;;

    return scale * R_pow(-log(unif_rand()), -1.0/shape);
}

double minvweibull(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	shape <= 0.0 ||
	order >= shape)
	return R_NaN;;

    return R_pow(scale, order) * gammafn(1.0 - order / shape);
}

double levinvweibull(double limit, double shape, double scale, double order, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	shape <= 0.0 ||
	order >= shape ||
	limit <= 0.0)
	return R_NaN;;

    tmp1 = R_pow(scale / limit, shape);
    tmp2 = order / shape;

    return R_pow(scale, order) * gammafn(1.0 - tmp2) * pgamma(tmp1, 1.0 - tmp2, 1.0, 0, 0) + R_pow(limit, order) * (1.0 - exp(-tmp1));
}
