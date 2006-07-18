/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse gamma distribution, to calculate raw moments and limited moments 
 *  of the random variable and to simulate random variates. See ../R/invgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgamma(double x, double shape, double scale, int give_log)
{
	double tmp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale)  ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

     if (!R_FINITE(x) || x < 0.0) 
      return R_D_d0;

	tmp = 1.0 / x;
    
    return  give_log ?
      -2.0 * log(x) + dgamma(tmp, shape, 1.0 / scale, 1) :
      R_pow(x, -2.0) * dgamma(tmp, shape, 1.0 / scale, 0);
}

double pinvgamma(double q, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;
  
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = 1.0 / q;
    
    return (lower_tail ? R_D_exp(pgamma(tmp, shape, 1.0 / scale, 0,1)):
	    R_D_exp(pgamma(tmp, shape, 1.0 / scale, 1,1)));
}

double qinvgamma(double p, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);

    return (lower_tail ? 1.0 / qgamma(tmp, shape, 1.0 / scale, 0, 0) :
	    1.0 / qgamma(tmp, shape, 1.0 / scale, 1, 1));
}

double rinvgamma(double shape, double scale)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    /*    return 1.0 / qgamma(unif_rand(), shape, 1.0 / scale, 1, 0); */
    return 1.0 / rgamma(shape, 1.0 / scale);
}

double minvgamma(double k, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(k) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	k >= shape)
	error(_("invalid arguments"));

    return R_pow(scale, k) * gammafn(shape - k) / gammafn(shape);
}

double levinvgamma(double d, double shape, double scale, double order, int give_log)
{
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(d) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order >= shape ||
	d <= 0.0)
	error(_("invalid arguments"));

    return R_pow(scale, order) * gammafn(shape - order) * (pgamma(1.0 / d, shape - order, 1.0 / scale, 0, 0)) / gammafn(shape) + R_pow(d, order) * pgamma(1.0 / d, shape, 1.0 / scale, 1, 0);
}
