/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the transformed beta distribution, and to simulate random
 *  variates. See ../R/trbeta.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dtrgamma(double x, double shape1, double scale, double shape2, int give_log)
{ 
  double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    tmp = R_pow(x, shape2);
    
    return  (give_log ? dgamma(tmp, shape1, scale, 1) + log(shape2) + (shape2 - 1.0) * log(x) :
	     shape2 * R_pow(x, shape2 - 1.0) * dgamma(tmp, shape1, scale, 0));
}

double ptrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0)
	error(_("invalid arguments"));

    tmp = R_pow(q, shape2);

    if (q <= 0)
	return R_DT_0;
    
    return (lower_tail ? R_D_exp(pgamma(tmp, shape1, scale, 1, 1)):
	    R_D_exp(pbeta(tmp, shape1, scale, 0, 1)));
}

double qtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);

    return (lower_tail ? R_D_exp((1.0 / shape2) * qgamma(p, shape1, scale, 1, 1)) :
	    R_D_exp((1.0 / shape2) * qgamma(p, shape1, scale, 0, 1)));
}

double rtrgamma(double shape1, double scale, double shape2)
{
    double a;	

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	error(_("invalid arguments"));
	
    a = unif_rand();

    return R_pow(qgamma(a, shape1, scale, 1, 0), 1.0 / shape2);
}
