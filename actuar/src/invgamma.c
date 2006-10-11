/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Inverse Gamma distribution. See ../R/invgamma.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgamma(double x, double shape, double scale, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale)  ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp1 = 1.0 / x;
    tmp2 = 1.0 / scale;

    return give_log ?
	-2.0 * log(x) + dgamma(tmp1, shape, tmp2, 1) :
	R_pow(x, -2.0) * dgamma(tmp1, shape, tmp2, 0);
}

double pinvgamma(double q, double shape, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    tmp1 = 1.0 / q;
    tmp2 = 1.0 / scale;

    return lower_tail ?
	R_D_exp(pgamma(tmp1, shape, tmp2, 0, 1)):
	R_D_exp(pgamma(tmp1, shape, tmp2, 1, 1));
}

double qinvgamma(double p, double shape, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp1 = R_D_qIv(p);
    tmp2 = 1.0 / scale;

    return lower_tail ?
	1.0 / qgamma(tmp1, shape, tmp2, 0, 0) :
	1.0 / qgamma(tmp1, shape, tmp2, 1, 0);
}

double rinvgamma(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    return 1.0 / rgamma(shape, 1.0 / scale);
}

double minvgamma(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order >= shape)
	return R_NaN;;

    return R_pow(scale, order) * gammafn(shape - order) / gammafn(shape);
}

double levinvgamma(double limit, double shape, double scale, double order, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order >= shape ||
	limit <= 0.0)
	return R_NaN;;

    tmp1 = 1.0 / limit;
    tmp2 = 1.0 / scale;

    return R_pow(scale, order) * gammafn(shape - order) * pgamma(tmp1, shape - order, tmp2, 0, 0) / gammafn(shape) + R_pow(limit, order) * pgamma(tmp1, shape, tmp2, 1, 0);
}
