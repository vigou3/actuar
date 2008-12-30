/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Fonctions to calculate raw and limited moments for the lognormal
 *  distribution. See ../R/LognormalMoments.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mlnorm(double order, double logmean, double logsd, int give_log)
{
    if (!R_FINITE(logmean) ||
        !R_FINITE(logsd)   ||
        !R_FINITE(order)   ||
        logsd <= 0.0)
        return R_NaN;

    return exp(order * (logmean + 0.5 * order * R_pow_di(logsd, 2)));
}

double levlnorm(double limit, double logmean, double logsd, double order,
                int give_log)
{
    double u;

    if (!R_FINITE(logmean) ||
        !R_FINITE(logsd)   ||
        !R_FINITE(order)   ||
        logsd <= 0.0)
        return R_NaN;

    if (limit <= 0.0)
        return 0;

    u = (log(limit) - logmean)/logsd;

    return exp(order * (logmean + 0.5 * order * R_pow(logsd, 2.0))) *
        pnorm(u - order * logsd, 0., 1.0, 1, 0) +
        R_VG__0(limit, order) * pnorm(u, 0., 1.0, 0, 0);
}
