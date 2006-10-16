/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the loggamma
 *  distribution. See ../R/lgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dlgamma(double x, double shapelog, double ratelog, int give_log)
{
    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0)
	return R_NaN;;

    if (!R_FINITE(x) || x < 1.0)
	return R_D_d0;

    return R_D_exp(dgamma(log(x), shapelog, 1.0 / ratelog, 1) - log(x));
}

double plgamma(double q, double shapelog, double ratelog, int lower_tail,
	       int log_p)
{
    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0)
	return R_NaN;;

    if (q <= 1.0)
	return R_DT_0;

    return pgamma(log(q), shapelog, 1.0 / ratelog, lower_tail, log_p);
}

double qlgamma(double p, double shapelog, double ratelog, int lower_tail,
	       int log_p)
{
    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0)
	return R_NaN;;

    R_Q_P01_boundaries(p, 1, R_PosInf);
    p = R_D_qIv(p);

    return exp(qgamma(p, shapelog, 1.0 / ratelog, lower_tail, 0));
}

double rlgamma(double shapelog, double ratelog)
{
    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0)
	return R_NaN;;

    return exp(rgamma(shapelog, 1.0 / ratelog));
}

double mlgamma(double order, double shapelog, double ratelog, int give_log)
{
    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	!R_FINITE(order) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0 ||
	order >= ratelog)
	return R_NaN;

    return R_pow(1.0 - order / ratelog, -shapelog);
}

double levlgamma(double limit, double shapelog, double ratelog, double order,
		 int give_log)
{
    double u;

    if (!R_FINITE(shapelog) ||
	!R_FINITE(ratelog) ||
	!R_FINITE(limit) ||
	!R_FINITE(order) ||
	shapelog <= 0.0 ||
	ratelog <= 0.0 ||
	limit <= 0.0 ||
	order >= ratelog)
	return R_NaN;

    u = log(limit);

    return R_pow(1.0 - order / ratelog, -shapelog)
	* pgamma(u * (ratelog - order), shapelog, 1.0, 1, 0)
	+ R_VG__0(limit, order) * pgamma(u * ratelog, shapelog, 1.0, 0, 0);
}
