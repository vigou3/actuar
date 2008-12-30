/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse Weibull distribution. See ../R/InverseWeibull.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvweibull(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u * e^(-u) / x
     *
     *  with u = (scale/x)^shape.
     */

    double logu;

    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    /* handle also x == 0 here */
    if (!R_FINITE(x) || x <= 0.0)
        return R_D__0;

    logu = shape * (log(scale) - log(x));

    return R_D_exp(log(shape) + logu - exp(logu) - log(x));
}

double pinvweibull(double q, double shape, double scale, int lower_tail,
                   int log_p)
{
    double u;

    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    if (q <= 0)
        return R_DT_0;

    u = exp(shape * (log(scale) - log(q)));

    return R_DT_val(exp(-u));
}

double qinvweibull(double p, double shape, double scale, int lower_tail,
                   int log_p)
{
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale * R_pow(-log(R_D_Lval(p)), -1.0 / shape);
}

double rinvweibull(double shape, double scale)
{
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    return scale * R_pow(rexp(1.0), -1.0 / shape);
}

double minvweibull(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        !R_FINITE(order) ||
        scale <= 0.0 ||
        shape <= 0.0 ||
        order >= shape)
        return R_NaN;;

    return R_pow(scale, order) * gammafn(1.0 - order / shape);
}

double levinvweibull(double limit, double shape, double scale, double order,
                     int give_log)
{
    double u, tmp;

    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        !R_FINITE(order) ||
        scale <= 0.0 ||
        shape <= 0.0 ||
        order >= shape)
        return R_NaN;;

    if (limit <= 0.0)
        return 0;

    tmp = 1.0 - order / shape;

    u = exp(shape * (log(scale) - log(limit)));

    return R_pow(scale, order) * gammafn(tmp) * pgamma(u, tmp, 1.0, 0, 0)
        + R_VG__0(limit, order) * (0.5 - exp(-u) + 0.5);
}
