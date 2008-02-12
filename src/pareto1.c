/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the single-parameter Pareto distribution. See
 *  ../R/SingleParameterPareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dpareto1(double x, double shape, double min, int give_log)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < min)
        return R_D__0;

    return R_D_exp(log(shape) + shape * log(min) - (shape + 1.0) * log(x));
}

double ppareto1(double q, double shape, double min, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (q <= min)
        return R_DT_0;

    return R_DT_Cval(R_pow(min / q, shape));
}

double qpareto1(double p, double shape, double min, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    R_Q_P01_boundaries(p, min, R_PosInf);
    p = R_D_qIv(p);

    return min / R_pow(R_D_Cval(p), 1.0 / shape);
}

double rpareto1(double shape, double min)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    return min / R_pow(unif_rand(), 1.0 / shape);
}

double mpareto1(double order, double shape, double min, int give_log)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        min <= 0.0   ||
        order >= shape)
        return R_NaN;

    return shape * R_pow(min, order) / (shape - order);
}

double levpareto1(double limit, double shape, double min, double order,
                  int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    tmp = shape - order;

    return shape * R_pow(min, order) / tmp
        - order * R_pow(min, shape) / (tmp * R_pow(limit, tmp));
}
