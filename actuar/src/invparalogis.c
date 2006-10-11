/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse paralogistic distribution. See ../R/invparalogis.R
 *  for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvparalogis(double x, double shape, double scale, int give_log)
{
    double tmp1, tmp2, tmp3;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp1 = x / scale;
    tmp2 = R_pow(tmp1, shape);
    tmp3 = shape * shape;

    return  give_log ?
	2.0 * log(shape) + tmp3 * (log(x) - log(scale)) - log(x) - (shape + 1.0) * log(1.0 + tmp2) :
	tmp3 * tmp2 / (x * R_pow(1.0 + tmp2, shape + 1.0));

}

double pinvparalogis(double q, double shape, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
	return 1;

    tmp1 = R_pow(scale, shape);
    tmp2 = R_pow(q, shape);

    return lower_tail ?
	R_D_exp(shape * (shape * log(q) - log(tmp1 + tmp2))):
	R_D_exp(log(1.0 - exp(shape * (shape * log(q) - log(tmp1 + tmp2)))));
}

double qinvparalogis(double p, double shape, double scale, int lower_tail, int log_p)
{
    double tmp, tmp1;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);

    tmp1 = -1.0 / shape;

    return lower_tail ?
	scale * R_pow(R_pow(tmp, tmp1) - 1.0, tmp1) :
	scale * R_pow(R_pow(1.0 - tmp, tmp1) - 1.0, tmp1);
}

double rinvparalogis(double shape, double scale)
{
    double a, tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    tmp = -1.0 / shape;

    a = unif_rand();

    return scale * R_pow(R_pow(a, tmp) - 1.0, tmp);
}

double minvparalogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -R_pow(shape, 2.0) ||
	order >= shape)
	return R_NaN;;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(shape + tmp) * gammafn(1.0 - tmp) / gammafn(shape);
}

double levinvparalogis(double limit, double shape, double scale, double order, int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -R_pow(shape, 2.0) ||
	order >= shape ||
	limit <= 0.0)
	return R_NaN;;

    tmp1 = R_pow(limit / scale, shape);
    tmp2 = order / shape;
    u = tmp1 / (1.0 + tmp1);

    return R_pow(scale, order) * gammafn(shape + tmp2) * gammafn(1.0 - tmp2) * pbeta(u, shape + tmp2, 1.0 - tmp2, 1, 0) + R_pow(limit, order) * (1.0 - R_pow(u, shape));
}
