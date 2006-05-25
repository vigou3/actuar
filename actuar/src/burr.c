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

double dburr(double x, double shape1, double scale, double shape2, int give_log)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
/*   return  (give_log ? 
 * 	     dbeta(x, 1.0, shape1, 1) + log(shape2) - log(scale) + (shape2 - 1.0)*(log(x) - log(scale)) + 2.0 * (-log(1.0 + exp(shape2 * (log(x) - log(scale))) *)) : 
 * 	     dbeta(x, 1.0, shape1, 0) * (shape2 / scale) * R_pow(x / scale, shape2 - 1.0) * R_pow(1.0 / (1.0 + R_pow(x / scale), shape2), 2.0)); 
 */

    return  (give_log ?
	log(shape1) + log(shape2) + shape2 * (log(x) - log(scale)) - log(x) - (shape1 + 1.0) * (log(R_pow(scale, shape2) + R_pow(x, shape2)) - shape2 * log(scale)) :
	shape1 * shape2 * R_pow(x / scale, shape2) / (x * R_pow(1.0 + R_pow(x / scale, shape2), shape1 + 1.0)));

}

double pburr(double x, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    if (x <= 0)
	return R_DT_0;

    return (lower_tail ? R_D_exp(log(1.0 - exp(shape1 * (shape2 * log(scale) - log(R_pow(scale, shape2) + R_pow(x, shape2)))))):
	    R_D_exp(shape1 * (shape2 * log(scale) - log(R_pow(scale, shape2) + R_pow(x, shape2)))));
}

double qburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
        shape2 <= 0.0)
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);

    return (lower_tail ? R_D_exp(log(scale) + (1.0 / shape2) * log(R_pow(1.0 / (1.0 - p), 1.0 / shape2) - 1.0)) :
	    R_D_exp(log(scale) + (1.0 / shape2) * log(R_pow(1.0 / p, 1.0 / shape2) - 1.0)));
}

double rburr(double shape1, double scale, double shape2)
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
	
    return scale * R_pow((1.0 - R_pow(a, (1.0 / shape1))) / (R_pow(a, (1.0 / shape1)) ), 1.0 / shape2);
}