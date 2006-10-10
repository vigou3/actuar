/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse Burr distribution. See ../R/invburr.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvburr(double x, double shape1, double shape2, double scale, int give_log)
{
    double tmp1, tmp2, tmp3;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp1 = x / scale;
    tmp2 = R_pow(tmp1, shape2);
    tmp3 = shape2 * (log(x) - log(scale));

    return give_log ?
	log(shape1) + log(shape2) + shape1 * tmp3 - log(x) - (shape1 + 1.0) * log(1.0 + tmp2) :
	shape1 * shape2 * R_pow(tmp2, shape1) / (x * R_pow(1 + tmp2, shape1 + 1.0));
}

double pinvburr(double q, double shape1, double shape2, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2, tmp3;

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

    tmp1 = q / scale;
    tmp2 = R_pow(tmp1, shape2);
    tmp3 = shape2 * (log(q) - log(scale));

    return lower_tail ?
	R_D_exp(shape1 * (tmp3 - log(1.0 + tmp2))):
	R_D_exp(log(1.0 - exp(shape1 * (tmp3 - log(1.0 + tmp2)))));
}

double qinvburr(double p, double shape1, double shape2, double scale, int lower_tail, int log_p)
{
    double tmp, tmp1;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);

    tmp1 = R_pow(tmp, 1.0 / shape1);

    return lower_tail ? scale * R_pow(tmp1 / (1.0 - tmp1), 1.0 / shape2) :
	scale * R_pow((1.0 - tmp1) / tmp1, 1.0 / shape2);
}

double rinvburr(double shape1, double shape2, double scale)
{
    double a, tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;

    a = unif_rand();
    tmp = R_pow(a, 1.0 / shape1);

    return scale * R_pow(tmp / (1.0 - tmp), 1.0 / shape2);
}

double minvburr(double order, double shape1, double shape2, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	order <= - shape1 * shape2 ||
	order >= shape2)
	return R_NaN;

    tmp = order / shape2;

    return R_pow(scale, order) * gammafn(shape1 + tmp) * gammafn(1.0 - tmp) / gammafn(shape1);
}

double levinvburr(double limit, double shape1, double shape2, double scale, double order, int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	limit <= 0.0 ||
	order <= -shape1 * shape2 ||
	order >= shape2)
	return R_NaN;

    tmp1 = limit / scale;
    tmp2 = order / shape2;
    u = R_pow(tmp1, shape2) / (1.0 + R_pow(tmp1, shape2));

    return R_pow(scale, order) * gammafn(shape1 + tmp2) * gammafn(1.0 - tmp2) * pbeta(u, shape1 + tmp2, 1.0 - tmp2, 1, 0) / gammafn(shape1) + R_pow(limit, order) * (1.0 - R_pow(u, shape1));
}
