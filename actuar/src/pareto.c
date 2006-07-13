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

    if (!R_FINITE(x) || x < 0.0) 
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
  double tmp, tmp1;

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

double mpareto(double k, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	!R_FINITE(k) ||
	shape <= 0.0 || 
	scale <= 0.0 ||
	k <= -1.0 ||
	k >= shape)
	error(_("invalid arguments"));

    return R_pow(scale, k) * gammafn(k + 1.0) * gammafn(shape - k) / gammafn(shape);

}

double levpareto(double x, double shape, double scale, double order, int give_log)
{
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	!R_FINITE(x) ||
	!R_FINITE(order) ||
	shape <= 0.0 || 
	scale <= 0.0 ||
	x <= 0.0)
	error(_("invalid arguments"));

    return R_pow(scale, order) * gammafn(order + 1.0) * gammafn(shape - order) * pbeta(x / (x + scale), order + 1.0, shape - order, 1, 0) / gammafn(shape) + R_pow(x, order) * R_pow(scale / (scale + x), shape);
}
