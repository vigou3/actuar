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

double dgenpareto(double x, double shape1, double scale, double shape2, int give_log)
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
      -log(beta(shape1, shape2)) + shape1 * log(scale) + (shape2 - 1.0) * log(x) - (shape1 + shape2) * log(x + scale) :
      (1 / beta(shape1, shape2)) * R_pow(scale, shape1) * R_pow(x, shape2 - 1.0) / R_pow(x + scale, shape1 + shape2);
}

double rgenpareto(double shape1, double scale, double shape2)
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

    return scale * ((qbeta(a, shape2, shape1, 1, 0)) / (1.0 - qbeta(a, shape2, shape1, 1, 0)));
}
