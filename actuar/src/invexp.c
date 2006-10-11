/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse exponential distribution. See ../R/invexp.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvexp(double x, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp = scale / x;

    return  give_log ?
	log(scale) - tmp - 2.0 * log (x) :
	tmp * exp(-tmp) / x;
}

double pinvexp(double q, double scale, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    tmp = scale / q;

    return lower_tail ? R_D_exp(-tmp): R_D_exp(log(1.0 - exp(-tmp)));
}

double qinvexp(double p, double scale, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);

    return lower_tail ? -scale / log(tmp) : -scale / log(1.0 - tmp);
}


double rinvexp(double scale)
{
    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    return -scale / log(unif_rand());
}

double minvexp(double order, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	order >= 1.0)
	return R_NaN;

    return R_pow(scale, order) * gammafn(1.0 - order);
}

double levinvexp(double limit, double scale, double order, int give_log)
{
    double tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	R_FINITE(order) ||
	scale <= 0.0 ||
	order <= 0.0 ||
	order >= 1.0 ||
	limit <= 0.0)
	return R_NaN;

    tmp = scale / limit;

    return R_pow(scale, order) * gammafn(1.0 - order) * pgamma(tmp, 1.0 - order, 1.0, 0, 0) + R_pow(limit, order) * (1.0 - exp(-tmp));
}
