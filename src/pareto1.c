/*  ===== actuar: An R Package for Actuarial Science =====
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
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(min))
	return x + shape + min;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < min)
        return ACT_D__0;

    return ACT_D_exp(log(shape) + shape * log(min) - (shape + 1.0) * log(x));
}

double ppareto1(double q, double shape, double min, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape) || ISNAN(min))
	return q + shape + min;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (q <= min)
        return ACT_DT_0;

    return ACT_DT_Cval(R_pow(min / q, shape));
}

double qpareto1(double p, double shape, double min, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape) || ISNAN(min))
	return p + shape + min;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, min, R_PosInf);
    p = ACT_D_qIv(p);

    return min / R_pow(ACT_D_Cval(p), 1.0 / shape);
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
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape) || ISNAN(min))
	return order + shape + min;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (order >= shape)
	return R_PosInf;

    return shape * R_pow(min, order) / (shape - order);
}

double levpareto1(double limit, double shape, double min, double order,
                  int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape) || ISNAN(min) || ISNAN(order))
	return limit + shape + min + order;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(min)   ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        min <= 0.0)
        return R_NaN;

    if (limit <= min)
        return 0.0;

    double tmp = shape - order;

    return shape * R_pow(min, order) / tmp
        - order * R_pow(min, shape) / (tmp * R_pow(limit, tmp));
}
