/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions, raw and limited moments and to simulate random variates
 *  for the inverse Pareto distribution. See ../R/invpareto.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvpareto(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u^shape * (1 - u) / x
     *
     *  with u = x/(x + scale) = 1/(1 + scale/x).
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    logu = - log(1 + scale / x);
    log1mu = - log(1 + x / scale);

    return R_D_exp(log(shape) + shape * logu + log1mu - log(x));
}

double pinvpareto(double q, double shape, double scale, int lower_tail,
		  int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    u = exp(-log1p(exp(log(scale) - log(q))));

    return R_DT_val(R_pow(q, shape));
}

double qinvpareto(double p, double shape, double scale, int lower_tail,
		  int log_p)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale / (R_pow(R_D_Lval(p), -1.0 / shape) - 1.0);
}

double rinvpareto(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    return scale / (R_pow(unif_rand(), -1.0 / shape) - 1.0);
}

double minvpareto(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= 1.0)
	return R_NaN;;

    return R_pow(scale, order) * gammafn(shape + order) * gammafn(1.0 - order)
	/ gammafn(shape);
}

double levinvpareto(double limit, double shape, double scale, double order,
		    int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape)
	return R_NaN;;

    if (limit <= 0.0)
	return 0;

    tmp1 = shape + order;
    tmp2 = 1.0 - order;

    u = exp(-log1p(exp(log(scale) - log(limit))));

    return R_pow(scale, order) * gammafn(tmp1) * gammafn(tmp2)
	* pbeta(u, tmp1, tmp2, 1, 0) / gammafn(shape)
	+ R_pow(limit, order) * (0.5 - R_pow(u, shape) + 0.5);
}
