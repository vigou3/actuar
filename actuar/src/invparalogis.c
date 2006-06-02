/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse paralogistic distribution, and to simulate random
 *  variates. See ../R/invparalogis.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvparalogis(double x, double shape, double scale, int give_log)
{

  double tmp1;
  double tmp2;
  double tmp3;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    tmp1 = R_pow(shape, 2.0);
    tmp2 = log(x) - log(scale);
    tmp3 = x / scale;
    
    return  give_log ?
	2.0 * log(shape) + tmp1 * tmp2 - log(x) - (shape + 1.0) * log(1.0 + exp(shape * tmp2)) :
  tmp1 * R_pow(tmp3, tmp1) / (x * R_pow(1.0 + R_pow(tmp3, shape), shape + 1.0));

}

double pinvparalogis(double q, double shape, double scale, int lower_tail, int log_p)
{
  double tmp1;
  double tmp2;
  
    if (!R_FINITE(shape) || 
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
        return 1;

	tmp1 = R_pow(shape, 2.0);
	tmp2 = log(q) - log(scale);
    
    return (lower_tail ? R_D_exp(tmp1 * tmp2 - shape * log(1.0 + exp(shape * tmp2))):
	    R_D_exp(log(1.0 - exp(tmp1 * tmp2 - shape * log(1.0 + exp(shape * tmp2))))));
}

double qinvparalogis(double p, double shape, double scale, int lower_tail, int log_p)
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

  tmp1 = 1 / shape;

  return (lower_tail ? scale * R_pow(R_pow(1.0 / (1.0 - tmp), tmp1) - 1.0, tmp1):
	    scale * R_pow(R_pow(1.0 / tmp, tmp1) - 1.0, tmp1));
}


double rinvparalogis(double shape, double scale)
{
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow((R_pow(a, 1.0 / shape)) / (1.0 - R_pow(a, 1.0 / shape)), 1.0 / shape);
}
