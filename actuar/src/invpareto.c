/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse Pareto distribution, to calculate raw moments of the random
 *  variable and to simulate random variates. See ../R/invpareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvpareto(double x, double shape, double scale, int give_log)
{
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      shape <= 0.0 || 
      scale <= 0.0) 
    return R_NaN;
  
  if (!R_FINITE(x) || x < 0.0) 
    return R_D_d0;
  
  return  give_log ?
    log(shape) + log (scale) + (shape - 1.0) * log(x) - (shape + 1.0) * log(x + scale) :
    shape * scale * R_pow(x, shape - 1.0) / R_pow(x + scale, shape + 1.0);  
}

double pinvpareto(double q, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;
  
  if (!R_FINITE(shape) || 
      !R_FINITE(scale) ||
      shape <= 0.0 || 
      scale <= 0.0)
    return R_NaN;;
  
  if (q <= 0)
    return R_DT_0;
  
  if (!R_FINITE(q))
    return 1;
  
  tmp = shape * (log(q) - log(q + scale));
  
  return (lower_tail ? R_D_exp(tmp):
	  R_D_exp(log(1.0 - exp(tmp))));
}

double qinvpareto(double p, double shape, double scale, int lower_tail, int log_p)
{
  double tmp;
  
  if (!R_FINITE(shape) || 
      !R_FINITE(scale) ||
      shape <= 0.0 || 
      scale <= 0.0)
    return R_NaN;;
  
  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);
  
  return (lower_tail ? scale * R_pow(tmp, 1.0 / shape) / (1.0 - R_pow(tmp, 1.0 / shape)) :
	  scale * R_pow(1.0 - tmp, 1.0 / shape) / (1.0 - R_pow(1.0 - tmp, 1.0 / shape))) ;
}

double rinvpareto(double shape, double scale)
{
  double a;
  
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      shape <= 0.0 ||
      scale <= 0.0)
    return R_NaN;;
  
  a = unif_rand();
  
  return scale * R_pow(a, 1.0 / shape) / (1.0 - R_pow(a, 1.0 / shape));
}

double minvpareto(double k, double shape, double scale, int give_log)
{
  
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      !R_FINITE(k) ||
      shape <= 0.0 ||
      scale <= 0.0 ||
      k <= -shape ||
      k >= 1.0)
    return R_NaN;;
  
  return R_pow(scale, k) * gammafn(shape + k) * gammafn(1.0 - k) / gammafn(shape);
}

double levinvpareto(double d, double shape, double scale, double order, int give_log)
{
  if (!R_FINITE(shape) || 
      !R_FINITE(scale) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      shape <= 0.0 || 
      scale <= 0.0 ||
      d <= 0.0 ||
      order <= -shape ||
      order >= 1.0)
    return R_NaN;;
  
  return R_pow(scale, order) * shape * gammafn(order + shape) * gammafn(1.0 - order) * pbeta(d / (d + scale), order + shape, 1.0 - order, 1, 0) / gammafn(shape + 1.0) + R_pow(d, order) * (1.0 - R_pow(d / (scale + d), shape));
}
