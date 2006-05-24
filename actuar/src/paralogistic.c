/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the Pareto distribution, and to simulate random
 *  variates. See ../R/pareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"

double dparalogistic(double x, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
      2.0 * log(shape) + 2.0 * (log(x) - log(scale)) - log(x) - (shape + 1.0) * log(1.0 + R_pow(x / scale, shape)) :
      R_pow(shape, 2.0) * R_pow(x / scale, shape) / (x * R_pow(1.0 + R_pow(x / scale, shape), shape + 1.0));
}

double rparalogistic(double shape, double scale)
{
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a =  unif_rand();
	
    return scale * R_pow(a / (1.0 - a), 1.0 / shape);
}
