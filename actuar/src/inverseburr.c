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

double dinverseburr(double x, double tau, double scale, double gamma, int give_log)
{
    if (!R_FINITE(tau) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	tau <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(tau) + log(gamma) + gamma * tau * (log(x) - log(scale)) - log(X) - (tau + 1.0)*log(1.0 + R_pow(1.0 + x / scale, gamma)) :
	tau * gamma * R_pow(x / scale, gamma * tau) / (x * R_pow(1 + R_pow(x / scale, gamma), tau + 1.0));
    
}

double rinverseburr(double tau, double scale, double gamma)
{	
    double a;

    if (!R_FINITE(tau) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	tau <= 0.0 ||
	scale <= 0.0 ||
	gamma <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow(R_pow(a, 1.0 / tau) / (1.0 - R_pow(a, 1.0 / tau)), 1.0 / gamma);
}
