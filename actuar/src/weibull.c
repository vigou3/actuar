/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to calculate raw moments and limited moments 
 *  of the Weibull distribution. See ../R/weibull.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mweibull(double k, double scale, double shape, int give_log)
{
  
  if (!R_FINITE(scale) ||
      !R_FINITE(shape) ||
      !R_FINITE(k) ||
      scale <= 0.0 ||
      shape <= 0.0 ||
      k <= -shape)
    return R_NaN;
  
  return R_pow(scale, k) * gammafn(1.0 + k / shape);
}

double levweibull(double d, double scale, double shape, double order, int give_log)
{
  
  if (!R_FINITE(scale) ||
      !R_FINITE(shape) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      scale <= 0.0 ||
      shape <= 0.0 ||
      d <= 0.0 ||
      order <= -shape)
    return R_NaN;
  
  return R_pow(scale, order) * gammafn(1.0 + order / shape) * pgamma(R_pow(d, shape), 1.0 + order / shape, 1.0 / scale, 1, 0) + R_pow(d, order) * exp(-R_pow(d / scale, shape));
}
