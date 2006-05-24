</*  ===== actuar: an R package for Actuarial Science =====
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

double dloglogistic(double x, double gamma, double scale, int give_log)
{
    if (!R_FINITE(gamma) ||
	!R_FINITE(scale) ||
	gamma <= 0.0 || 
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(gamma) + gamma * (log(x) - log(scale)) - log(x) - 2.0 * log(1 + R_pow(x / scale, gamma)) :
	gamma * R_pow(x / scale, gamma)/(x * R_pow(1.0 + R_pow(x / scale, gamma), 2.0));
}

double rloglogistic(double gamma, double scale)
{	
    double a;
	
    if (!R_FINITE(gamma) ||
	!R_FINITE(scale) ||
	gamma <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow(a / (1.0 - a), 1.0 / gamma);
}
