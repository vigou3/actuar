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
    /*  We work with the density expressed as
     *
     *  u * e^(-u) / x
     *
     *  with u = scale/x.
     */

    double logu;

    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    logu = log(scale) - log(x);

    return R_D_exp(logu - exp(logu) - log(x));
}

double pinvexp(double q, double scale, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    u = exp(log(scale) - log(q));

    return R_DT_val(exp(-u));
}

double qinvexp(double p, double scale, int lower_tail, int log_p)
{
    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return -scale / log(R_D_Lval(p));
}

double rinvexp(double scale)
{
    if (!R_FINITE(scale) || scale <= 0.0)
	return R_NaN;

    return -scale / rexp(1.0);
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
    double u, tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	R_FINITE(order) ||
	scale <= 0.0 ||
	order <= 0.0 ||
	order >= 1.0 ||
	limit <= 0.0)
	return R_NaN;

    tmp = 1.0 - order;

    u = exp(log(scale) - log(limit));

    return R_pow(scale, order) * gammafn(tmp) * pgamma(u, tmp, 1.0, 0, 0)
	+ R_pow(limit, order) * (0.5 - exp(-u) + 0.5);
}
