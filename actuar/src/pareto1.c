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

double dpareto1(double x, double shape, double min, int give_log)
{
    if (!R_FINITE(shape) ||
	(!R_FINITE(min) ||
	shape <= 0.0 || 
	min <= 0.0 ||
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(shape) + shape * log(min) - (shape + 1.0) * log(x):
	shape * R_pow(min, shape) / R_pow(x, shape + 1.0)

double rpareto1(double shape, double min)
{	
    if (!R_FINITE(shape) ||
	!R_FINITE(min) ||
	shape <= 0.0 ||
	min <= 0.0)
	error(_("invalid arguments"));

    return min / R_pow(unif_rand(), 1.0 / shape);
}
