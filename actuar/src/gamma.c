/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Gamma
 *  distribution.
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
	scale <= 0.0 ||
	order <= -shape)
	return R_NaN;

    return R_pow(scale, order) * gammafn(order + shape) / gammafn(shape);
}

double levgamma(double limit, double shape, double scale, double order,
		int give_log)
{
    double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	limit <= 0.0 ||
	order <= -shape)
	return R_NaN;

    tmp = order + shape;

    return R_pow(scale, order) * gammafn(tmp) *
	pgamma(limit, tmp, scale, 1, 0) / gammafn(shape) +
	R_VG__0(limit, order) * pgamma(limit, shape, scale, 0, 0);
}
