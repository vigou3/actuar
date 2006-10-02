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

double dinvtrgamma(double x, double shape1, double scale, double shape2, int give_log)
{
  double tmp1, tmp2;
  
  if (!R_FINITE(shape1) ||
      !R_FINITE(scale)  ||
      !R_FINITE(shape2) ||
      shape1 <= 0.0 || 
      scale <= 0.0 ||
      shape2 <= 0.0) 
    return R_NaN;
  
  if (!R_FINITE(x) || x < 0.0) 
    return R_D_d0;
  
  tmp1 = R_pow(1.0 / x, shape2);
  tmp2 = R_pow(scale, -shape2);

  return  give_log ?
    log(shape2) - (shape2 + 1.0) * log(x) + dgamma(tmp1, shape1, tmp2, 1) :
    shape2 *  R_pow(x, -shape2 - 1.0) * dgamma(tmp1, shape1, tmp2, 0);
}

double pinvtrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
  double tmp1, tmp2;
  
  if (!R_FINITE(shape1) ||
      !R_FINITE(scale) ||
      !R_FINITE(shape2) ||
      shape1 <= 0.0 || 
      scale <= 0.0 ||
      shape2 <= 0.0) 
    return R_NaN;;
  
  if (q <= 0)
    return R_DT_0;
  
  tmp = R_pow(1.0 / q, shape2);
  tmp2 = R_pow(scale, -shape2);
  
  return (lower_tail ? R_D_exp(pgamma(tmp1, shape1, tmp2, 0,1)):
	  R_D_exp(pgamma(tmp1, shape1, tmp2, 1,1)));
}

double qinvtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
  double tmp1, tmp2, tmp3;
  
  if (!R_FINITE(shape1) ||
      !R_FINITE(scale) ||
      !R_FINITE(shape2) ||
      shape1 <= 0.0 || 
      scale <= 0.0 ||
      shape2 <= 0.0)
    return R_NaN;;
  
  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp1 = R_D_qIv(p);
  tmp2 = R_pow(scale, -shape2);
  tmp3 = 1.0 / shape2
  
  return (lower_tail ? R_pow(qgamma(tmp1, shape1, tmp2, 0, 0), -tmp3) :
	  R_pow(qgamma(tmp1, shape1, tmp2, 1, 0), -tmp3));
}

double rinvtrgamma(double shape1, double scale, double shape2)
{
  double a;
  
  if (!R_FINITE(shape1) ||
      !R_FINITE(scale) ||
      !R_FINITE(shape2) ||
      shape1 <= 0.0 || 
      scale <= 0.0 ||
      shape2 <= 0.0)
    return R_NaN;;
  
  a = rgamma(shape1, 1.0);
  
  return R_pow(a, -1.0 / shape2) / scale;
}

double minvtrgamma(double k, double shape1, double scale, double shape2, int give_log)
{
  
  if(!R_FINITE(shape1) ||
     !R_FINITE(scale) ||
     !R_FINITE(shape2) ||
     !R_FINITE(k) ||
     shape1 <= 0.0 || 
     scale <= 0.0 ||
     shape2 <= 0.0 ||
     k >= shape1 * shape2)
    return R_NaN;;
  
  return R_pow(scale, k) * gammafn(shape1 - k / shape2) / gammafn(shape1);
}

double levinvtrgamma(double d, double shape1, double scale, double shape2, double order, int gve_log)
{
  double u, tmp2;
  
  if (!R_FINITE(shape1) ||
      !R_FINITE(scale) ||
      !R_FINITE(shape2) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      shape1 <= 0.0 || 
      scale <= 0.0 ||
      shape2 <= 0.0 ||
      order >= shape1 * shape2 ||
      d <= 0.0)
    return R_NaN;;
  
  u = R_pow(1.0 / d, shape2);
  tmp2 = R_pow(scale, -shape2)
  
  return R_pow(scale, order) * gammafn(shape1 - order / shape2) * pgamma(u, shape1 - order / shape2, tmp2, 0, 0) / gammafn(shape1) + R_pow(d, order) * (pgamma(u, shape1, tmp2, 1, 0)) ;
}
