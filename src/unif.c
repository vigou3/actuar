/*  ===== actuar: An R Package for Actuarial Science =====
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
    double tmp;

    if (!R_FINITE(min) ||
        !R_FINITE(max) ||
        min >= max)
        return R_NaN;

    if (order == -1.0)
        return (log(fabs(max)) - log(fabs(min))) / (max - min);

    tmp = order + 1;

    return (R_pow(max, tmp) - R_pow(min, tmp)) / ((max - min) * tmp);
}

double levunif(double limit, double min, double max, double order, int give_log)
{
    double tmp;

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

    tmp = order + 1;

    return (R_pow(limit, tmp) - R_pow(min, tmp)) / ((max - min) * tmp) +
        R_pow(limit, order) * (max - limit) / (max - min);
}

double mgfunif(double x, double min, double max, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(min) ||
        !R_FINITE(max) ||
        min >= max)
        return R_NaN;

    if (x == 0.0)
        return ACT_D_exp(0.0);

    tmp1 = exp(x * max) - exp(x * min);
    tmp2 = x * (max - min);

    return ACT_D_exp(log(tmp1) - log(tmp2));
}
