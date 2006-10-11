/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse transformed gamma distribution. See
 *  ../R/invtrgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvtrgamma(double x, double shape1, double shape2, double scale, int give_log)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale)  ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    tmp1 = R_pow(x, -shape2);
    tmp2 = R_pow(scale, -shape2);

    return  give_log ?
	log(shape2) - (shape2 + 1.0) * log(x) + dgamma(tmp1, shape1, tmp2, 1) :
	shape2 *  tmp1 * dgamma(tmp1, shape1, tmp2, 0) / x;
}

double pinvtrgamma(double q, double shape1, double shape2, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    tmp1 = R_pow(q, -shape2);
    tmp2 = R_pow(scale, -shape2);

    return lower_tail ?
	R_D_exp(pgamma(tmp1, shape1, tmp2, 0, 1)):
	R_D_exp(pgamma(tmp1, shape1, tmp2, 1, 1));
}

double qinvtrgamma(double p, double shape1, double shape2, double scale, int lower_tail, int log_p)
{
    double tmp1, tmp2, tmp3;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp1 = R_D_qIv(p);
    tmp2 = R_pow(scale, -shape2);
    tmp3 = -1.0 / shape2;

    return lower_tail ?
	R_pow(qgamma(tmp1, shape1, tmp2, 0, 0), tmp3) :
	R_pow(qgamma(tmp1, shape1, tmp2, 1, 0), tmp3);
}

double rinvtrgamma(double shape1, double shape2, double scale)
{
    double a;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	return R_NaN;;

    a = rgamma(shape1, 1.0);

    return scale * R_pow(a, -1.0 / shape2);
}

double minvtrgamma(double order, double shape1, double shape2, double scale, int give_log)
{

    if(!R_FINITE(shape1) ||
       !R_FINITE(scale) ||
       !R_FINITE(shape2) ||
       !R_FINITE(order) ||
       shape1 <= 0.0 ||
       scale <= 0.0 ||
       shape2 <= 0.0 ||
       order >= shape1 * shape2)
	return R_NaN;;

    return R_pow(scale, order) * gammafn(shape1 - order / shape2) / gammafn(shape1);
}

double levinvtrgamma(double limit, double shape1, double shape2, double scale, double order, int give_log)
{
    double u, tmp1, tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	order >= shape1 * shape2 ||
	limit <= 0.0)
	return R_NaN;;

    u = R_pow(limit, -shape2);
    tmp1 = R_pow(scale, -shape2);
    tmp2 = order / shape2;

    return R_pow(scale, order) * gammafn(shape1 - tmp2) * pgamma(u, shape1 - tmp2, tmp1, 0, 0) / gammafn(shape1) + R_pow(limit, order) * pgamma(u, shape1, tmp1, 1, 0);
}
