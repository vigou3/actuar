/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to calculate raw moments and limited moments 
 *  of the Gamma distribution. See ../R/gamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mgamma(double k, double shape, double scale, int give_log)
{
  
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      !R_FINITE(k) ||
      shape <= 0.0 ||
      scale <= 0.0 ||
      k <= -shape)
    return R_NaN;
  
  return R_pow(scale, k) * gammafn(k + shape) / gammafn(shape);
}

double levgamma(double d, double shape, double scale, double order, int give_log)
{
  
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      shape <= 0.0 ||
      scale <= 0.0 ||
      d <= 0.0 ||
      order <= -shape)
    return R_NaN;
  
  return R_pow(scale, order) * gammafn(shape + order) * pgamma(d, shape + order, 1.0 / scale, 1, 0) / gammafn(shape) + R_pow(d, order) * pgamma(d, shape, 1.0 / scale, 0, 0);
}
