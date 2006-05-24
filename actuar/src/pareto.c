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

double dpareto(double x, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) || !R_FINITE(scale) ||
	shape <= 0.0 || scale <= 0.0 || x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(shape) + shape * log(scale) - (shape + 1.0) * log(x + scale) :
	shape * R_pow(scale, shape) / R_pow(x + scale, shape + 1.0);
}

double rpareto(double shape, double scale)
{
    if (!R_FINITE(shape) || !R_FINITE(scale) ||
	shape <= 0.0 || scale <= 0.0)
	error(_("invalid arguments"));

    return scale * (R_pow(unif_rand(), -1.0 / shape) - 1.0);
}
