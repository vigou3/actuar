/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Fonctions to calculate raw and limited moments for the Weibull
 *  distribution. See ../R/WeibullMoments.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mweibull(double order, double shape, double scale, int give_log)
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

    if (order <= -shape)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(1.0 + order / shape);
}

double levweibull(double limit, double shape, double scale, double order,
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

    if (order <= -shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u, tmp;

    tmp = 1.0 + order / shape;

    u = exp(shape * (log(limit) - log(scale)));

    return R_pow(scale, order) * gammafn(tmp) * pgamma(u, tmp, 1.0, 1, 0) +
        ACT_DLIM__0(limit, order) * exp(-u);
}
