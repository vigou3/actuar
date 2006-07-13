/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the loglogistic distribution, and to simulate random
 *  variates. See ../R/llogis.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dllogis(double x, double shape, double scale, int give_log)
{
  double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

     if (!R_FINITE(x) || x < 0.0) 
      return R_D_d0;

    tmp1 = shape * (log(x) - log(scale));
    tmp2 = R_pow(x / scale, shape);
    
    return  give_log ?
	log(shape) + tmp1 - log(x) - 2.0 * log(1.0 + exp(tmp1)) :
	shape * tmp2 / (x * R_pow(1.0 + tmp2, 2.0));
}

double pllogis(double q, double shape, double scale, int lower_tail, int log_p)
{

  double tmp1, tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
        return 1;

    tmp1 = shape + (log(q) - log(scale));
    tmp2 = R_pow(q / scale, shape);

    return (lower_tail ? R_D_exp(shape * log(q) - log(R_pow(scale, shape) + R_pow(q, shape))):
	    R_D_exp(log(1.0 - exp(shape * log(q) - log(R_pow(scale, shape) + R_pow(q, shape))))));
}

double qllogis(double p, double shape, double scale, int lower_tail, int log_p)
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

    return (lower_tail ? scale * R_pow(tmp / (1.0 - tmp), tmp1) :
	    scale * R_pow((1.0 - tmp) / tmp, tmp1));
}

double rllogis(double shape, double scale)
{	
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow(a / (1.0 - a), 1.0 / shape);
}

double mllogis(double k, double shape, double scale, int give_log)
{	
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(k) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	k <= -shape ||
	k >= shape)
	error(_("invalid arguments"));

    return R_pow(scale, k) * gammafn(1.0 + k / shape) * gammafn(1.0 - k / shape);
}

double levllogis(double x, double shape, double scale, double order, int give_log)
{	
  double temp1, temp2;
       
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	!R_FINITE(x) ||
	!R_FINITE(order) ||
	shape <= 0.0 ||
	scale <= 0.0 ||
	order <= -shape ||
	x <= 0.0)
	error(_("invalid arguments"));

    temp1 = R_pow(x / scale, shape);
    temp2 = temp1 / (1.0 + temp1);

    return R_pow(scale, order) * gammafn(1.0 + order / shape) * gammafn(1.0 - order / shape) * pbeta(temp2, 1.0 + order / shape, 1.0 - order / shape, 1, 0) + R_pow(x, order) * (1.0 - temp2);
}
