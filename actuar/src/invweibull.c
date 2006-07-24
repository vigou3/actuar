/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse Weibull distribution, to calculate raw moments and limited moments 
 *  of the random variable and to simulate random variates. See ../R/invweibull.R for details.
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
  
  if (!R_FINITE(scale) ||
      !R_FINITE(shape)  ||
      scale <= 0.0 || 
      shape <= 0.0) 
    return R_NaN;;
  
  if (!R_FINITE(x) || x < 0.0) 
    return R_D_d0;
  
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
    return R_NaN;;
  
  if (q <= 0)
    return R_DT_0;
  
  tmp = R_pow(scale / q, shape);
  
  return (lower_tail ? R_D_exp(-tmp):
	  R_D_exp(log(1.0 - exp(tmp))));
}

double qinvweibull(double p, double scale, double shape, int lower_tail, int log_p)
{
  double tmp, tmp1;
  
  if (!R_FINITE(scale) || 
      !R_FINITE(shape) ||
      scale <= 0.0 || 
      shape <= 0.0)
    return R_NaN;;
  
  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);
  
  tmp1 = 1 / shape;
  
  return (lower_tail ? scale * R_pow(log(1.0 / tmp), -tmp1):
	  scale * R_pow(log(1.0 / (1.0 - tmp)), -tmp1));
}

double rinvweibull(double scale, double shape)
{
  if (!R_FINITE(scale) ||
      !R_FINITE(shape) ||
      scale <= 0.0 ||
      shape <= 0.0)
    return R_NaN;;
  
  return shape * scale / log(1.0 / unif_rand());
}

double minvweibull(double k, double scale, double shape, int give_log)
{
  if (!R_FINITE(scale) ||
      !R_FINITE(shape) ||
      !R_FINITE(k) ||
      scale <= 0.0 ||
      shape <= 0.0 ||
      k >= shape)
    return R_NaN;;
  
  return R_pow(scale, k) * gammafn(1.0 - k / shape);
}

double levinvweibull(double d, double scale, double shape, double order, int give_log)
{
  if (!R_FINITE(scale) ||
      !R_FINITE(shape) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      scale <= 0.0 ||
      shape <= 0.0 ||
      order >= shape ||
      d <= 0.0)
    return R_NaN;;
  
  return R_pow(scale, order) * gammafn(1.0 - order / shape) * (1.0 - pgamma(R_pow(1.0 / d, shape), 1.0 - order / shape, 1.0 / scale, 1, 0)) + R_pow(d, order) * (1.0 - exp(-R_pow(scale / d, shape)));
}
