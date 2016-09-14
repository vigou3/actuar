/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Exponential
 *  distribution. See ../R/ExponentialSupp.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mexp(double order, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(scale))
	return order + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(order) ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -1.0)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(1.0 + order);
}

double levexp(double limit, double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(scale) || ISNAN(order))
	return limit + scale + order;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(order) ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -1.0)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u, tmp;

    tmp = 1.0 + order;
    u = exp(log(limit) - log(scale));

    return R_pow(scale, order) * gammafn(tmp) *
        pgamma(u, tmp, 1.0, 1, 0) +
        ACT_DLIM__0(limit, order) * exp(-u);
}

double mgfexp(double x, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(scale))
	return x + scale;
#endif
    if (!R_FINITE(scale) ||
        scale <= 0.0 ||
        scale * x > 1.0)
        return R_NaN;

    if (x == 0.0)
        return ACT_D_exp(0.0);

    return ACT_D_exp(-log1p(-scale * x));
}
