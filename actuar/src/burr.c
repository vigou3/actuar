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

double dburr(double x, double shape, double scale, double gamma, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  (give_log ?
	log(shape) + log(gamma) + gamma * (log(x) - log(scale)) - log(x) - (shape + 1.0) * (log(R_pow(scale, gamma) + R_pow(x, gamma)) - gamma * log(scale)) :
	shape * gamma * R_pow(x / scale, gamma) / (x * R_pow(1.0 + R_pow(x / scale, gamma), shape + 1.0)));
}

double pburr(double x, double shape, double scale, double gamma, int lower_tail, int log_p)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    if (lower_tail)
	return log_p ? 
	  
    else:  
        return log_p ?
	  
}

double rburr(double shape, double scale, double gamma)
{
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	gamma <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();
	
    return scale * R_pow((1.0 - R_pow(a, (1.0 / shape))) / (R_pow(a, (1.0 / shape)) ), 1.0 / gamma);
}
