/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Inverse Gamma distribution. See ../R/invgamma.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgamma(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  u^shape * e^(-u) / (x * gamma(shape))
     *
     *  with u = scale/x.
     */

    double logu;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return R_D_d0;

    logu = log(scale) - log(x);

    return R_D_exp(shape * logu - exp(logu) - log(x) - lgamma(shape));
}

double pinvgamma(double q, double shape, double scale, int lower_tail,
		 int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    if (q <= 0)
	return R_DT_0;

    u = exp(log(scale) - log(q));

    return pgamma(u, shape, 1.0, !lower_tail, log_p);
}

double qinvgamma(double p, double shape, double scale, int lower_tail,
		 int log_p)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 0, R_PosInf);
    p = R_D_qIv(p);

    return scale / qgamma(p, shape, 1.0, !lower_tail, 0);
}

double rinvgamma(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	return R_NaN;;

    return scale / rgamma(shape, 1.0);
}

double minvgamma(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order >= shape)
	return R_NaN;;

    return R_pow(scale, order) * gammafn(shape - order) / gammafn(shape);
}

double levinvgamma(double limit, double shape, double scale, double order,
		   int give_log)
{
    double u, tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order >= shape)
	return R_NaN;;

    if (limit <= 0.0)
	return 0;

    tmp = shape - order;

    u = exp(log(scale) - log(limit));

    return R_pow(scale, order) * gammafn(shape - order)
	* pgamma(u, tmp, 1.0, 0, 0) / gammafn(shape)
	+ R_VG__0(limit, order) * pgamma(u, shape, 1.0, 1, 0);
}
