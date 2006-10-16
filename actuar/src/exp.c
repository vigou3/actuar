/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Exponential
 *  distribution.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mexp(double order, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	order <= -1.0)
	return R_NaN;

    return R_pow(scale, order) * gammafn(1.0 + order);
}

double levexp(double limit, double scale, double order, int give_log)
{
    double tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	limit <= 0.0 ||
	order <= -1.0)
	return R_NaN;

    tmp = 1.0 + order;

    return R_pow(scale, order) * gammafn(tmp) *
	pgamma(limit, tmp, scale, 1, 0) +
	R_VG__0(limit, order) * exp(-limit * scale);
}
