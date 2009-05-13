/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the paralogistic distribution. See ../R/Paralogistic.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dparalogis(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape^2 * u^shape * (1 - u) / x
     *
     *  with u = 1/(1 + v), v = (x/scale)^shape.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return R_D__0;

    /* handle x == 0 separately */
    if (x == 0)
    {
	if (shape < 1) return R_PosInf;
	if (shape > 1) return R_D__0;
	/* else */
	return R_D_val(1.0 / scale);
    }

    tmp = shape * (log(x) - log(scale));
    logu = - log1p(exp(tmp));
    log1mu = - log1p(exp(-tmp));

    return R_D_exp(2.0 * log(shape) + shape * logu + log1mu - log(x));
}

double pparalogis(double q, double shape, double scale, int lower_tail,
                  int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (q <= 0)
        return R_DT_0;

    u = exp(-log1p(exp(shape * (log(q) - log(scale)))));

    return R_DT_Cval(R_pow(u, shape));
}

double qparalogis(double p, double shape, double scale, int lower_tail,
                  int log_p)
{
    double tmp;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    tmp = 1.0 / shape;

    return scale * R_pow(R_pow(R_D_Cval(p), -tmp) - 1.0, tmp);
}

double rparalogis(double shape, double scale)
{
    double tmp;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    tmp = 1.0 / shape;

    return scale * R_pow(R_pow(unif_rand(), -tmp) - 1.0, tmp);
}

double mparalogis(double order, double shape, double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
	return R_NaN;

    if (order <= -shape ||
        order >= shape * shape)
        return R_PosInf;

    tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(shape - tmp)
        / gammafn(shape);
}

double levparalogis(double limit, double shape, double scale, double order,
                    int give_log)
{
    double u, tmp1, tmp2, tmp3;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    tmp1 = order / shape;
    tmp2 = 1.0 + tmp1;
    tmp3 = shape - tmp1;

    u = exp(-log1p(exp(shape * (log(limit) - log(scale)))));

    return R_pow(scale, order) * gammafn(tmp2) * gammafn(tmp3)
        * pbeta(0.5 - u + 0.5, tmp2, tmp3, 1, 0) / gammafn(shape)
        + R_VG__0(limit, order) * R_pow(u, shape);
}
