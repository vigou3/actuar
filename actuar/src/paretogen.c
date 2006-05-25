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

double dparetogen(double x, double shape, double scale, double tau, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(tau) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	tau <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
      -log(beta(shape, tau)) + shape * log(scale) + (tau - 1.0) * log(x) - (shape + tau) * log(x + scale) :
      (1 / beta(shape, tau)) * R_pow(scale, shape) * R_pow(x, tau - 1.0) / R_pow(x + scale, shape + tau);
}

double rparetogen(double shape, double scale, double tau)
{	
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(tau) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	tau <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * ((qbeta(a, tau, shape, 1, 0)) / (1.0 - qbeta(a, tau, shape, 1, 0)));
}
