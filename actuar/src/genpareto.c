/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Generalized Pareto distribution.. See ../R/genpareto.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dgenpareto(double x, double shape1, double shape2, double scale,
		  int give_log)
{
    /*  We work with the density expressed as
     *
     *  scale^shape1 * u^shape2 * (1 - u)^shape1 / (x * beta(shape1, shape2)
     *
     *  with u = x/(x + scale).
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

    tmp = x + scale;
    log1mu = - log(exp(tmp));
    logu = log(x) + log1mu;

    return R_D_exp(shape1 * log(scale) + shape2 * logu + shape1 * log1mu
		   - log(x) - lbeta(shape2, shape1));
}

double pgenpareto(double q, double shape1, double shape2, double scale,
		  int lower_tail, int log_p)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
	return 1;

    return pbeta(q / (q + scale), shape2, shape1, lower_tail, log_p);
}

double qgenpareto(double p, double shape1, double shape2, double scale,
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

    return scale / (1.0 - 1.0 / qbeta(p, shape2, shape1, lower_tail, 0));
}

double rgenpareto(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    return scale / (1.0 - 1.0 / rbeta(shape2, shape1));
}

double mgenpareto(double order, double shape1, double shape2, double scale,
		  int give_log)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	!R_FINITE(order)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= -shape2 ||
	order  >= shape1)
	return R_NaN;

    return R_pow(scale, order) * beta(shape1 - order, shape2 + order)
	/ beta(shape1, shape2);
}

double levgenpareto(double limit, double shape1, double shape2, double scale,
		    double order, int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(scale)  ||
	!R_FINITE(order)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= -shape2)
	return R_NaN;

    if (limit <= 0.0)
	return 0;

    if (!R_FINITE(limit))
	return mgenpareto(order, shape1, shape2, scale, 0);

    tmp1 = shape1 - order;
    tmp2 = shape2 + order;

    u = limit / (limit + scale);

    return R_pow(scale, order) * beta(tmp1, tmp2) / beta(shape1, shape2)
	* pbeta(u, tmp2, tmp1, 1, 0)
	+ R_pow(limit, order) * pbeta(u, shape2, shape1, 0, 0);
}
