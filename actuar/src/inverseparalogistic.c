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

double dinverseparalogistic(double x, double tau, double scale, int give_log)
{
    if (!R_FINITE(tau) ||
	!R_FINITE(scale) ||
	tau <= 0.0 || 
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	2.0 * log(tau) + R_pow(tau, 2.0) * (log(x) - log(scale)) - log(x) - (tau + 1.0) * log(1.0 + R_pow(x / scale, tau))) :
    R_pow(tau, 2.0) * R_pow(x / scale, R_pow(tau, 2.0)) / (x * R_pow(1.0 + R_pow(x / scale, tau), tau + 1.0));

    
}

double rinverseparalogistic(double tau, double scale)
{
    double a;
	
    if (!R_FINITE(tau) ||
	!R_FINITE(scale) ||
	tau <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow((R_pow(a, 1.0 / tau)) / (1.0 - R_pow(a, 1.0 / tau)), 1.0 / tau);
}
