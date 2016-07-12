/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Pareto distribution. See ../R/Pareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dpareto(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u^shape * (1 - u) / x
     *
     *  with u = 1/(1 + v), v = x/scale.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0) return ACT_D_val(shape / scale);

    tmp = log(x) - log(scale);
    logu = - log1pexp(tmp);
    log1mu = - log1pexp(-tmp);

    return ACT_D_exp(log(shape) + shape * logu + log1mu - log(x));
}

double ppareto(double q, double shape, double scale, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    u = exp(-log1pexp(log(q) - log(scale)));

    return ACT_DT_Cval(R_pow(u, shape));
}

double qpareto(double p, double shape, double scale, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * (R_pow(ACT_D_Cval(p), -1.0 / shape) - 1.0);
}

double rpareto(double shape, double scale)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * (R_pow(unif_rand(), -1.0 / shape) - 1.0);
}

double mpareto(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0  ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -1.0 ||
        order >= shape)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(1.0 + order) * gammafn(shape - order)
        / gammafn(shape);
}

double levpareto(double limit, double shape, double scale, double order,
                 int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -1.0)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    tmp1 = 1.0 + order;
    tmp2 = shape - order;

    u = exp(-log1pexp(log(limit) - log(scale)));

    return R_pow(scale, order) * gammafn(tmp1) * gammafn(tmp2)
        * pbeta(0.5 - u + 0.5, tmp1, tmp2, 1, 0) / gammafn(shape)
        + ACT_DLIM__0(limit, order) * R_pow(u, shape);
}
