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

double dinvweibull(double x, double scale, double shape, int give_log)
{
  
  double tmp;

    if (R_FINITE(scale) ||
	R_FINITE(shape)  ||
	scale <= 0.0 || 
	shape <= 0.0 ||
	x < 0.0) 
	error(_("invalid arguments"));

    tmp = R_pow(scale / x, shape);
    
    return  give_log ?
      log(shape) + shape * (log(scale) - log(x)) - tmp - log(x) :
      shape * tmp * exp(-tmp) / x;
}

double pinvweibull(double q, double scale, double shape, int lower_tail, int log_p)
{
  double tmp;
  
    if (!R_FINITE(scale) || 
	!R_FINITE(shape) ||
	scale <= 0.0 || 
	shape <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = R_pow(scale / q, shape);
    
    return (lower_tail ? R_D_exp(-tmp):
	    R_D_exp(log(1.0 - exp(tmp))));
}

double qinvweibull(double p, double scale, double shape, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(scale) || 
	!R_FINITE(shape) ||
	scale <= 0.0 || 
	shape <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, 1);

  tmp = 1 / shape;

    return (lower_tail ? R_D_exp(log(scale) - tmp * log(log(1 / p))):
	    R_D_exp(log(scale) - tmp  * log(log(1 / (1 - p)))));
}

double rinvweibull(double scale, double shape)
{
    if (!R_FINITE(scale) ||
	!R_FINITE(shape) ||
	scale <= 0.0 ||
	shape <= 0.0)
	error(_("invalid arguments"));

    return shape * scale / log(1.0 / unif_rand());
}
