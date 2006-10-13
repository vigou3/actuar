/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed beta distribution. See ../R/trbeta.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

/*  To take advantage of optimizations in R's dbeta.c, we code the
 *  density as
 *
 *      shape2 * u * (1 - u) * dbeta(u, shape3, shape1) / x,
 *
 *  with u = v/(1 + v), v = (x/scale)^shape2. Furthermore, in the
 *  hope to reduce rounding error and/or over/underflows, we work
 *  in the log of the terms.
 */

double dtrbeta(double x, double shape1, double shape2, double shape3,
	       double scale, int give_log)
{
    double logv, logu, log1mu, tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale) ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    logv = shape2 * (log(x) - log(scale));
    log1mu = - log1p(exp(logv));
    logu = logv + log1mu;

    tmp = dbeta(exp(logu), shape3, shape1, 1) +
	logu + log1mu + log(shape2) - log(x);

    return give_log ? tmp : exp(tmp);
}

double ptrbeta(double q, double shape1, double shape2, double shape3,
	       double scale, int lower_tail, int log_p)
{
    double tmp, u;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
	return 1;

    tmp = shape2 * (log(q) - log(scale));
    u = exp(tmp - log1p(exp(tmp)));

    return pbeta(u, shape3, shape1, lower_tail, log_p);
}

double qtrbeta(double p, double shape1, double shape2, double shape3,
	       double scale, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale  <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    tmp = qbeta(p, shape3, shape1, lower_tail, 0);

    return scale * R_pow(tmp / (1.0 - tmp), 1.0 / shape2);
}

double rtrbeta(double shape1, double shape2, double shape3, double scale)
{
    double a;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale) ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    a = rbeta(shape3, shape1);

    return scale * R_pow(a / (1.0 - a), 1.0 / shape2);
}

double mtrbeta(double order, double shape1, double shape2, double shape3,
	       double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale)  ||
	!R_FINITE(order)  ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= - shape3 * shape2 ||
	order  >= shape1 * shape2)
	return R_NaN;

    tmp = order / shape2;

    return R_pow(scale, order) * beta(shape3 + tmp, shape1 - tmp) /
	beta(shape1, shape3);
}

double levtrbeta(double limit, double shape1, double shape2, double shape3,
		 double scale, double order, int give_log)
{
    double u, tmp, tmp1, tmp2, tmp3;

    if (!R_FINITE(shape1) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0 ||
	scale  <= 0.0 ||
	order  <= - shape3 * shape2)
	return R_NaN;

    if (!R_FINITE(limit))
	return mtrbeta(order, shape1, shape2, shape3, scale, 0);

    tmp1 = order / shape2;
    tmp2 = shape3 + tmp1;
    tmp3 = shape1 - tmp1;

    tmp = shape2 * (log(limit) - log(scale));
    u = exp(tmp - log1p(exp(tmp)));

    return R_pow(scale, order) * beta(tmp2, tmp3) / beta(shape1, shape3) *
	pbeta(u, tmp2, tmp3, 1, 0) +
	R_pow(limit, order) * pbeta(u, shape3, shape1, 0, 0);
}
