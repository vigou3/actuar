/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to calculate raw and limited moments for the Weibull
 *  distribution. See ../R/weibull.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mweibull(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	shape <= 0.0 ||
	order <= -shape)
	return R_NaN;

    return R_pow(scale, order) * gammafn(1.0 + order / shape);
}

double levweibull(double d, double shape, double scale, double order, int give_log)
{
    double u, tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	shape <= 0.0 ||
	limit <= 0.0 ||
	order <= -shape)
	return R_NaN;

    u = R_pow(limit / scale, shape);
    tmp = 1.0 + order / shape;

    return R_pow(scale, order) * gammafn(tmp) * pgamma(u, tmp, 1.0, 1, 0) +
	R_pow(limit, order) * exp(-u);
}
