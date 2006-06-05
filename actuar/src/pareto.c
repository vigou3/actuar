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

double dpareto(double x, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

    if (!R_FINITE(x)  ||
	x < 0.0) 
      return R_D_d0;
    
    return (give_log ?
	log(shape) + shape * log(scale) - (shape + 1.0) * log(x + scale) :
	shape * R_pow(scale, shape) / R_pow(x + scale, shape + 1.0));
}

double ppareto(double q, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;
  
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;
    
    tmp = log(scale) - log(q + scale);
    
    return (lower_tail ? R_D_exp(log(1.0 - exp(shape * (tmp)))):
	    R_D_exp(shape * (tmp)));
}

double qpareto(double p, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;
  double tmp1;

  if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);
  tmp1 = 1.0 / shape;

    return (lower_tail ? scale * (R_pow(1.0 - tmp, -tmp1) - 1.0) :
	    scale * (R_pow(tmp, -tmp1) - 1.0));
}

double rpareto(double shape, double scale)
{
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    return scale * (R_pow(unif_rand(), -1.0 / shape) - 1.0);
}
