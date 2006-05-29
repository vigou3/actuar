/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse gamma distribution, and to simulate random
 *  variates. See ../R/invgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgamma(double x, double scale, double shape, int give_log)
{
	double tmp;

    if (R_FINITE(scale) ||
	R_FINITE(shape)  ||
	scale <= 0.0 || 
	shape <= 0.0 ||
	x < 0.0) 
	error(_("invalid arguments"));

	tmp = 1.0 / x;
    
    return  give_log ?
      -2.0 * log(x) + dgamma(tmp, shape, scale, 1) :
      R_pow(x, -2.0) * dgamma(tmp, shape, scale, 0);
}

double pinvgamma(double q, double scale, double shape, int lower_tail, int log_p)
{
  double tmp;
  
    if (!R_FINITE(scale) || 
	!R_FINITE(shape) ||
	scale <= 0.0 || 
	shape <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = 1.0 / q;
    
    return (lower_tail ? R_D_exp(pgamma(tmp, shape, scale, 0,1)):
	    R_D_exp(pgamma(tmp, shape, scale, 1,1)));
}

double qinvgamma(double p, double scale, double shape, int lower_tail, int log_p)
{

  if (!R_FINITE(scale) || 
	!R_FINITE(shape) ||
	scale <= 0.0 || 
	shape <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, 1);

    return (lower_tail ? R_D_exp(-qgamma(1.0 - p, shape, scale, 1, 1)) :
	    R_D_exp(-qgamma(1.0 - p, shape, scale, 0, 1)));
}

double rinvgamma(double scale, double shape)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	error(_("invalid arguments"));

    return qgamma(unif_rand(), shape, scale, 1, 0);
}
