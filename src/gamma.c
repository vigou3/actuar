/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to calculate raw and limited moments for the Gamma
 *  distribution. See ../R/GammaSupp.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mgamma(double order, double shape, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape) || ISNAN(scale))
	return order + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(order + shape) / gammafn(shape);
}

double levgamma(double limit, double shape, double scale, double order,
                int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape) || ISNAN(scale) || ISNAN(order))
	return limit + shape + scale + order;
#endif
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

    double u, tmp;

    tmp = order + shape;
    u = exp(log(limit) - log(scale));

    return R_pow(scale, order) * gammafn(tmp) *
        pgamma(u, tmp, 1.0, 1, 0) / gammafn(shape) +
        ACT_DLIM__0(limit, order) * pgamma(u, shape, 1.0, 0, 0);
}

double mgfgamma(double t, double shape, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(t) || ISNAN(shape) || ISNAN(scale))
	return t + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0 ||
        scale * t > 1.)
        return R_NaN;

    if (t == 0.0)
        return ACT_D__1;

    return ACT_D_exp(-shape * log1p(-scale * t));
}
