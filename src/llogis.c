/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the loglogistic distribution. See ../R/Loglogistic.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dllogis(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u * (1 - u) / x
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    tmp = shape * (log(x) - log(scale));
    logu = - log1p(exp(-tmp));
    log1mu = - log1p(exp(tmp));

    return R_D_exp(log(shape) + logu + log1mu - log(x));
}

double pllogis(double q, double shape, double scale, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (q <= 0)
	return R_DT_0;

    u = exp(-log1p(exp(shape * (log(scale) - log(q)))));

    return R_DT_val(u);
}

double qllogis(double p, double shape, double scale, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale * R_pow(1.0 / R_D_Cval(p) - 1.0, 1.0/shape);
}

double rllogis(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    return scale * R_pow(1.0 / unif_rand() - 1.0, 1.0 / shape);
}

double mllogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	order >= shape)
	return R_NaN;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(1.0 - tmp);
}

double levllogis(double limit, double shape, double scale, double order,
		 int give_log)
{
    double u, tmp1, tmp2, tmp3;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape)
	return R_NaN;

    tmp1 = order / shape;
    tmp2 = 1.0 + tmp1;
    tmp3 = 1.0 - tmp1;

    u = exp(-log1p(exp(shape * (log(scale) - log(limit)))));

    return R_pow(scale, order) * gammafn(tmp2) * gammafn(tmp3)
	* pbeta(u, tmp2, tmp3, 1, 0)
	+ R_VG__0(limit, order) * (0.5 - u + 0.5);
}
