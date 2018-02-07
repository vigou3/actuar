/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to calculate raw and limited moments for the Uniform
 *  distribution. See ../R/UniformSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double munif(double order, double min, double max, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(min) || ISNAN(max))
	return order + min + max;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(max) ||
        min >= max)
        return R_NaN;

    if (order == -1.0)
        return (log(fabs(max)) - log(fabs(min))) / (max - min);

    double tmp = order + 1;

    return (R_pow(max, tmp) - R_pow(min, tmp)) / ((max - min) * tmp);
}

double levunif(double limit, double min, double max, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(min) || ISNAN(max) || ISNAN(order))
	return limit + min + max + order;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(max) ||
        min >= max)
        return R_NaN;

    if (limit <= min)
        return R_pow(limit, order);

    if (limit >= max)
        return munif(order, min, max, give_log);

    if (order == -1.0)
        return (log(fabs(limit)) - log(fabs(min))) / (max - min) +
            (max - limit) / (limit * (max - min));

    double tmp = order + 1;

    return (R_pow(limit, tmp) - R_pow(min, tmp)) / ((max - min) * tmp) +
        R_pow(limit, order) * (max - limit) / (max - min);
}

double mgfunif(double t, double min, double max, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(t) || ISNAN(min) || ISNAN(max))
	return t + min + max;
#endif
    if (!R_FINITE(min) ||
        !R_FINITE(max) ||
        min >= max)
        return R_NaN;

    if (t == 0.0)
        return ACT_D__1;

    double tmp1, tmp2;

    tmp1 = exp(t * max) - exp(t * min);
    tmp2 = t * (max - min);

    return ACT_D_exp(log(tmp1) - log(tmp2));
}
