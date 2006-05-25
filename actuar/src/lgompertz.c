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

double dlgompertz(double x, double scale, double shape, int give_log)
{
    if (R_FINITE(scale) ||
	R_FINITE(shape)  ||
	scale <= 0.0 || 
	shape <= 0.0 ||
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
      log(shape) + shape * (log(scale) - log(x)) - R_pow(scale / x, shape) - log(x) :
      shape * R_pow(scale / x, shape) * exp(-R_pow(scale / x, shape)) / x;
}

double rlgompertz(double scale, double shape)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	error(_("invalid arguments"));

    return shape * scale / log(1.0 / unif_rand());
}
