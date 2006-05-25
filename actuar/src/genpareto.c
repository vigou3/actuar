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
      dbeta(x, shape2, shape1, 1) + log(scale) - 2.0 * log(x + scale) :
      dbeta(x, shape2, shape1, 0) * scale / R_pow(x + scale,2.0);
}

double pgenpareto(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    u = q / (q + scale);

    if (q <= 0)
	return R_DT_0;
    
    return (lower_tail ? R_D_exp(pbeta(u, shape2, shape1, 1, 0)):
	    R_D_exp(pbeta(u, shape2, shape1, 0, 0)));
}

double qgenpareto(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);

    return (lower_tail ? R_D_exp(log(scale) + qbeta(p, shape2, shape1, 1, 1) - log(1.0 - qbeta(p, shape2, shape1, 1, 0))) :
	    R_D_exp(log(scale) + qbeta(p, shape2, shape1, 0, 1) - log(1.0 - qbeta(p, shape2, shape1, 0, 0))));
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
