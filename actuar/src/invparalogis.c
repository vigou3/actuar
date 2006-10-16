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
    /*  We work with the density expressed as
     *
     *  shape^2 * u^shape * (1 - u) / x
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp = shape * (log(x) - log(scale));
    logu = - log1p(exp(-tmp));
    log1mu = - log1p(exp(tmp));

    return R_D_exp(2.0 * log(shape) + shape * logu + log1mu - log(x));
}

double pinvparalogis(double q, double shape, double scale, int lower_tail,
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

    u = exp(-log1p(exp(shape * (log(scale) - log(q)))));

    return R_DT_val(R_pow(u, shape));
}

double qinvparalogis(double p, double shape, double scale, int lower_tail,
		     int log_p)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    tmp = -1.0 / shape;

    return scale * R_pow(R_pow(R_D_Lval(p), tmp) - 1.0, tmp);
}

double rinvparalogis(double shape, double scale)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    tmp = -1.0 / shape;

    return scale * R_pow(R_pow(unif_rand(), tmp) - 1.0, tmp);
}

double minvparalogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= - shape * shape ||
	order >= shape)
	return R_NaN;;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(shape + tmp) * gammafn(1.0 - tmp)
	/ gammafn(shape);
}

double levinvparalogis(double limit, double shape, double scale, double order,
		       int give_log)
{
    double u, tmp1, tmp2, tmp3;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape * shape)
	return R_NaN;;

    tmp1 = order / shape;
    tmp2 = shape + tmp1;
    tmp3 = 1.0 - tmp1;

    u = exp(-log1p(exp(shape * (log(scale) - log(limit)))));

    return R_pow(scale, order) * gammafn(tmp2) * gammafn(tmp3)
	* pbeta(u, tmp2, tmp3, 1, 0) / gammafn(shape)
	+ R_VG__0(limit, order) * (0.5 - R_pow(u, shape) + 0.5);
}
