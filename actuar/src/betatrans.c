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
	     dbeta(x, tau, shape, 1) + log(gamma) - log(scale) + (gamma - 1.0)*(log(x) - log(scale)) + 2.0 * (-log(1.0 + exp(gamma * (log(x) - log(scale))))) :
	     dbeta(x, tau, shape, 0) * (gamma / scale) * R_pow(x / scale, gamma - 1.0) * R_pow(1.0 / (1.0 + R_pow(x / scale), gamma), 2.0));
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
	tau <= 0.0) 
	error(_("invalid arguments"));

    u = R_pow(x / scale, gamma) / (1.0 + R_pow(x / scale, gamma));

    if (x <= 0)
	return R_DT_0;
    
    return (lower_tail ? R_D_exp(pbeta(u, tau, shape, 1, 0)):
	    R_D_exp(pbeta(u, tau, shape, 0, 0)));
}

double qbetatrans(double p, double shape, double scale, double gamma, double tau, int lower_tail, int log_p)
{

  if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(gamma) ||
	!R_FINITE(tau) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	gamma <= 0.0 || 
	tau <= 0.0) 
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);

    return (lower_tail ? R_D_exp(log(scale) + gamma * (qbeta(p, tau, shape, 1, 1)) - gamma * log(1.0 - qbeta(p, tau, shape, 1, 0))) :
	    R_D_exp(log(scale) + gamma * (1beta(p, tau, shape, 0, 1)) - gamma * log(1.0 - qbeta(p, tau, shape, 0, 0))));
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