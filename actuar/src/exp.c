/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to calculate raw moments and limited moments 
 *  of the Exponential distribution. See ../R/exp.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mexp(double k, double rate, int give_log)
{
  
  if (!R_FINITE(rate) ||
      !R_FINITE(k) ||
      rate <= 0.0 ||
      k <= -1.0)
    return R_NaN;
  
  return R_pow(rate, k) * gammafn(1.0 + k);
}

double levexp(double d, double rate, double order, int give_log)
{
  
  if (!R_FINITE(rate) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      rate <= 0.0 ||
      d <= 0.0 ||
      order <= -1.0)
    return R_NaN;
  
  return R_pow(rate, order) * gammafn(1.0 + order) * pgamma(d, order + 1.0, 1.0 / rate, 1, 0) + R_pow(d, order) * exp(-d / rate);
}
