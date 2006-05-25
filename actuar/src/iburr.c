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
#include "dpq.h"

double diburr(double x, double shape1, double scale, double shape2, int give_log)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(shape1) + log(shape2) + shape2 * shape1 * (log(x) - log(scale)) - log(x) - (shape1 + 1.0)*log(1.0 + R_pow(1.0 + x / scale, shape2)) :
	shape1 * shape2 * R_pow(x / scale, shape2 * shape1) / (x * R_pow(1 + R_pow(x / scale, shape2), shape1 + 1.0));
    
}

double riburr(double shape1, double scale, double shape2)
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

    return scale * R_pow(R_pow(a, 1.0 / shape1) / (1.0 - R_pow(a, 1.0 / shape1)), 1.0 / shape2);
}
