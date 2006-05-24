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

double dbetatrans(double x, double shape, double scale, double gamma, double tau, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	!R_FINITE(tau) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	tau <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  (give_log ?
	log(scale) + (1.0 / gamma) * (dbeta(x, tau, shape, 1) - log(1.0 - dbeta(x, tau, alpha, 0))) :
	scale * R_pow(dbeta(x, tau, shape, 0)/(1.0 - dbeta(x, tau, shape, 0)), 1.0 / gamma));
}

double pbetatrans(double q, double shape, double scale, double gamma, double tau, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	!R_FINITE(tau) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	tau <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    u = R_pow(x / scale, gamma) / (1.0 + R_pow(x / scale, gamma));
    
    if (lower_tail)
	return (log_p ? 
	  pbeta(u, tau, shape, 1, 1): 
	  pbeta(u, tau, shape, 1, 0));
    else:  
        return (log_p ?
	  pbeta(u, tau, shape, 0, 1):
          pbeta(u, tau, shape, 0, 0));
}

double rbetatrans(double shape, double scale, double gamma, double tau)
{
    double a;	

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	!R_FINITE(tau) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	gamma <= 0.0 ||
	tau <= 0.0)
	error(_("invalid arguments"));
	
    a = unif_rand();

    return scale * R_pow((qbeta(a, tau, shape, 1, 0) / (1.0 - qbeta(a, tau, shape, 1, 0))), (1.0 / gamma));
}
