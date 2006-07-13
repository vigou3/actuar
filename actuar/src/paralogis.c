/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the paralogistic distribution, and to simulate random
 *  variates. See ../R/paralogis.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dparalogis(double x, double shape, double scale, int give_log)
{
  double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

     if (!R_FINITE(x) || x < 0.0) 
      return R_D_d0;

    tmp1 = log(x) - log(scale);
    tmp2 = R_pow(x / scale, shape);
    
    return  give_log ?
      2.0 * log(shape) + shape * tmp1 - log(x) - (shape + 1.0) * log(1.0 + exp(shape * tmp1)) :
      R_pow(shape, 2.0) * tmp2 / (x * R_pow(1.0 + tmp2, shape + 1.0));
}

double pparalogis(double x, double shape, double scale, int lower_tail, int log_p)
{

  double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

    if (x <= 0)
	return R_DT_0;

    tmp1 = R_pow(scale, shape);
    tmp2 = R_pow(x, shape);

    return (lower_tail ? R_D_exp(log(1.0 - exp(shape * ( shape * log(scale) - log(tmp1 + tmp2))))):
	    R_D_exp(shape * ( shape * log(scale) - log(tmp1 + tmp2))));
}

double qparalogis(double p, double shape, double scale, int lower_tail, int log_p)
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

    return (lower_tail ? scale * R_pow(R_pow(1.0 - tmp, -tmp1) - 1.0, tmp1) :
	    scale * R_pow(R_pow(tmp, -tmp1) - 1.0, tmp1));
}

double rparalogis(double shape, double scale)
{
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a =  R_pow(unif_rand(), 1.0 / shape);
	
    return scale * R_pow((1.0 - a) / a, 1.0 / shape);
}

double mparalogis(double k, double shape, double scale, int give_log)
{
 
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(k) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	k <= -shape ||
	k >= R_pow(shape, 2.0))
	error(_("invalid arguments"));
	
    return R_pow(scale, k) * gammafn(1.0 + k / shape) * gammafn(shape - k / shape) / gammafn(shape);
}

double levparalogis(double d, double shape, double scale, double order, int give_log)
{
  double temp;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(d) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	d <= 0.0)
	error(_("invalid arguments"));

    temp = 1.0 / (1.0 + R_pow(d / scale, shape));
	
    return R_pow(scale, order) * gammafn(1.0 + order / shape) * gammafn(shape - order / shape) * pbeta(1.0 - temp, 1.0 + order / shape, shape - order / shape, 1, 0) + R_pow(d, order) * R_pow(temp, shape);
}
