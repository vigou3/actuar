/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed gamma distribution. See ../R/TransformedGamma.R
 *  for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dtrgamma(double x, double shape1, double shape2, double scale,
                int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape2 * u^shape1 * e^(-u) / (x * gamma(shape1))
     *
     *  with u = (x/scale)^shape2.
     */

    double logu;

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
    if (x == 0) R_D_mode(shape1 * shape2 > 1);

    logu = shape2 * (log(x) - log(scale));

    return R_D_exp(log(shape2) + shape1 * logu - exp(logu)
                   - log(x) - lgamma(shape1));
}

double ptrgamma(double q, double shape1, double shape2, double scale,
                int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return R_DT_0;

    u = exp(shape2 * (log(q) - log(scale)));

    return pgamma(u, shape1, 1.0, lower_tail, log_p);
}

double qtrgamma(double p, double shape1, double shape2, double scale,
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

    return scale * R_pow(qgamma(p, shape1, 1.0, lower_tail, 0),
                         1.0 / shape2);
}

double rtrgamma(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    return scale * R_pow(rgamma(shape1, 1.0), 1.0 / shape2);
}

double mtrgamma(double order, double shape1, double shape2, double scale,
                int give_log)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale <= 0.0 ||
        order <= -shape1 * shape2)
        return R_NaN;

    return R_pow(scale, order) * gammafn(shape1 + order/shape2)
        / gammafn(shape1);
}

double levtrgamma(double limit, double shape1, double shape2, double scale,
                  double order, int give_log)
{
    double u, tmp;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0 ||
        order <= -shape1 * shape2)
        return R_NaN;

    if (limit <= 0.0)
        return 0;

    tmp = shape1 + order / shape2;

    u = exp(shape2 * (log(limit) - log(scale)));

    return R_pow(scale, order) * gammafn(tmp)
        * pgamma(u, tmp, 1.0, 1, 0) / gammafn(shape1)
        + R_VG__0(limit, order) * pgamma(u, shape1, 1.0, 0, 0);
}
