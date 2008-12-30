/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Generalized Pareto distribution.. See ../R/GeneralizedPareto.R
 *  for details.
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
     *  u^shape2 * (1 - u)^shape1 / (x * beta(shape1, shape2))
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = x/scale.
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
        return R_D__0;

    /* handle x == 0 separately */
    if (x == 0) R_D_mode(shape2 > 1);

    tmp = log(x) - log(scale);
    logu = - log1p(exp(-tmp));
    log1mu = - log1p(exp(tmp));

    return R_D_exp(shape2 * logu + shape1 * log1mu - log(x)
                   - lbeta(shape2, shape1));
}

double pgenpareto(double q, double shape1, double shape2, double scale,
                  int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape1) ||
        !R_FINITE(scale)  ||
        !R_FINITE(shape2) ||
        shape1 <= 0.0 ||
        scale <= 0.0 ||
        shape2 <= 0.0)
        return R_NaN;

    if (q <= 0)
        return R_DT_0;

    u = exp(-log1p(exp(log(scale) - log(q))));

    return pbeta(u, shape2, shape1, lower_tail, log_p);
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

    return scale / (1.0 / qbeta(p, shape2, shape1, lower_tail, 0) - 1.0);
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

    return scale / (1.0 / rbeta(shape2, shape1) - 1.0);
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

    tmp1 = shape1 - order;
    tmp2 = shape2 + order;

    u = exp(-log1p(exp(log(scale) - log(limit))));

    return R_pow(scale, order) * beta(tmp1, tmp2) / beta(shape1, shape2)
        * pbeta(u, tmp2, tmp1, 1, 0)
        + R_VG__0(limit, order) * pbeta(u, shape2, shape1, 0, 0);
}
