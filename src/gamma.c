/*  ===== actuar: An R Package for Actuarial Science =====
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
        R_VG__0(limit, order) * pgamma(u, shape, 1.0, 0, 0);
}

double mgfgamma(double x, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0 ||
        scale * x > 1.)
        return R_NaN;

    if (x == 0.0)
        return R_D_exp(0.0);

    return R_D_exp(-shape * log1p(-scale * x));
}
