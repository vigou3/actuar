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

double dpareto1(double x, double shape, double min, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(min) ||
	shape <= 0.0 || 
	min <= 0.0)
	error(_("invalid arguments"));

    if (x <= min) error(_("invalid x"));
    
    return  give_log ?
	log(shape) + shape * log(min) - (shape + 1.0) * log(x):
      shape * R_pow(min, shape) / R_pow(x, shape + 1.0);
}

double ppareto1(double q, double shape, double min, int lower_tail, int log_p)
{
  double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(min) ||
	shape <= 0.0 ||
	min <= 0.0)
	error(_("invalid arguments"));

    if (q <= min)
	return R_DT_0;

    tmp = log(min) - log(q);
    
    return (lower_tail ? R_D_exp(log(1 - exp(shape * (tmp)))):
	    R_D_exp(shape * (tmp)));
}

double qpareto1(double p, double shape, double min, int lower_tail, int log_p)
{

  double tmp;

    R_Q_P01_boundaries(p, 0, 1);

    if (!R_FINITE(shape) ||
	!R_FINITE(min) ||
	shape <= 0.0 ||
	min <= 0.0)
	error(_("invalid arguments"));

	  tmp = (1.0 / shape);

    return (lower_tail ? R_D_exp(log(min) - tmp * log(1.0 - p)) :
	    R_D_exp(log(min) - tmp * log(p)));
}


double rpareto1(double shape, double min)
{	
    if (!R_FINITE(shape) ||
	!R_FINITE(min) ||
	shape <= 0.0 ||
	min <= 0.0)
	error(_("invalid arguments"));

    return min / R_pow(unif_rand(), 1.0 / shape);
}
