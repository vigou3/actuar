/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Exponential
 *  distribution. See ../R/ExponentialSupp.R for details.
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
    double u, tmp;

    if (!R_FINITE(scale) ||
	!R_FINITE(order) ||
	scale <= 0.0 ||
	order <= -1.0)
	return R_NaN;

    if (limit <= 0.0)
	return 0;

    tmp = 1.0 + order;

    u = exp(log(limit) - log(scale));

    return R_pow(scale, order) * gammafn(tmp) *
	pgamma(u, tmp, 1.0, 1, 0) +
	R_VG__0(limit, order) * exp(-u);
}

double mgfexp(double t, double scale, int give_log)
{
	/*check arguments */
	if (!R_FINITE(scale) ||
	    scale <= 0.0 ||
	    scale * t > 1.)
	  return R_NaN;
	  
	
	if(t == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( -log(1. - scale * t) ) ;
}