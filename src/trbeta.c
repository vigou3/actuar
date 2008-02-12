/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed beta distribution. See ../R/TransformedBeta.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dtrbeta(double x, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape2 * u^shape3 * (1 - u)^shape1 / (x * beta(shape1, shape3))
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape2.
     */

    double tmp, logu, log1mu;

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
        return R_D__0;

    tmp = shape2 * (log(x) - log(scale));
    logu = - log1p(exp(-tmp));
    log1mu = - log1p(exp(tmp));

    return R_D_exp(log(shape2) + shape3 * logu + shape1 * log1mu
                   - log(x) - lbeta(shape3, shape1));
}

double ptrbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
    double u;

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

    u = exp(-log1p(exp(-shape2 * (log(q) - log(scale)))));

    return pbeta(u, shape3, shape1, lower_tail, log_p);
}

double qtrbeta(double p, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
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

    return scale * R_pow(1.0 / qbeta(p, shape3, shape1, lower_tail, 0) - 1.0,
                         -1.0 / shape2);
}

double rtrbeta(double shape1, double shape2, double shape3, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * R_pow(1.0 / rbeta(shape3, shape1) - 1.0, -1.0 / shape2);
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

    return R_pow(scale, order) * beta(shape3 + tmp, shape1 - tmp)
        / beta(shape1, shape3);
}

double levtrbeta(double limit, double shape1, double shape2, double shape3,
                 double scale, double order, int give_log)
{
    double u, tmp1, tmp2, tmp3;

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

    if (limit <= 0.0)
        return 0;

    tmp1 = order / shape2;
    tmp2 = shape3 + tmp1;
    tmp3 = shape1 - tmp1;

    u = exp(-log1p(exp(-shape2 * (log(limit) - log(scale)))));

    return R_pow(scale, order) * beta(tmp2, tmp3) / beta(shape1, shape3)
        * pbeta(u, tmp2, tmp3, 1, 0)
        + R_VG__0(limit, order) * pbeta(u, shape3, shape1, 0, 0);
}
