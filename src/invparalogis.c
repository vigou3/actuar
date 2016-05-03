/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse paralogistic distribution. See ../R/InverseParalogistic.R
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
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0)
    {
	if (shape < 1.0) return R_PosInf;
	if (shape > 1.0) return ACT_D__0;
	/* else */
	return ACT_D_val(1.0 / scale);
    }

    tmp = shape * (log(x) - log(scale));
    logu = - log1pexp(-tmp);
    log1mu = - log1pexp(tmp);

    return ACT_D_exp(2.0 * log(shape) + shape * logu + log1mu - log(x));
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
        return ACT_DT_0;

    u = exp(-log1pexp(shape * (log(scale) - log(q))));

    return ACT_DT_val(R_pow(u, shape));
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

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    tmp = -1.0 / shape;

    return scale * R_pow(R_pow(ACT_D_Lval(p), tmp) - 1.0, tmp);
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
        scale <= 0.0)
        return R_NaN;

    if (order <= - shape * shape ||
        order >= shape)
	return R_PosInf;

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
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape * shape)
	return R_PosInf;

    tmp1 = order / shape;
    tmp2 = shape + tmp1;
    tmp3 = 1.0 - tmp1;

    u = exp(-log1pexp(shape * (log(scale) - log(limit))));

    return R_pow(scale, order) * gammafn(tmp2) * gammafn(tmp3)
        * pbeta(u, tmp2, tmp3, 1, 0) / gammafn(shape)
        + ACT_DLIM__0(limit, order) * (0.5 - R_pow(u, shape) + 0.5);
}
