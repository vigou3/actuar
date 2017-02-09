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
#include "actuar.h"

double dinvweibull(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u * e^(-u) / x
     *
     *  with u = (scale/x)^shape.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
	return x + shape + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    /* handle also x == 0 here */
    if (!R_FINITE(x) || x <= 0.0)
        return ACT_D__0;

    double logu = shape * (log(scale) - log(x));

    return ACT_D_exp(log(shape) + logu - exp(logu) - log(x));
}

double pinvweibull(double q, double shape, double scale, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape) || ISNAN(scale))
	return q + shape + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(shape * (log(scale) - log(q)));

    return ACT_DT_val(exp(-u));
}

double qinvweibull(double p, double shape, double scale, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape) || ISNAN(scale))
	return p + shape + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * R_pow(-log(ACT_D_Lval(p)), -1.0 / shape);
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
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape) || ISNAN(scale))
	return order + shape + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        !R_FINITE(order) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;

    if (order >= shape)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(1.0 - order / shape);
}

double levinvweibull(double limit, double shape, double scale, double order,
                     int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape) || ISNAN(scale) || ISNAN(order))
	return limit + shape + scale + order;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(shape) ||
        !R_FINITE(order) ||
        scale <= 0.0 ||
        shape <= 0.0)
        return R_NaN;

    if (order >= shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u = exp(shape * (log(scale) - log(limit)));

    return R_pow(scale, order) * actuar_gamma_inc(1.0 - order/shape, u)
        + ACT_DLIM__0(limit, order) * (0.5 - exp(-u) + 0.5);
}
