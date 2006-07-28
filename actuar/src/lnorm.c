/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to calculate raw moments and limited moments 
 *  of the Normal distribution. See ../R/normal.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mlnorm(double k, double mu, double sigma, int give_log)
{
  
  if (!R_FINITE(mu) ||
      !R_FINITE(sigma) ||
      !R_FINITE(k) ||
      sigma <= 0.0)
    return R_NaN;
  
  return exp(k * mu + 0.5 * R_pow(k, 2.0) * R_pow(sigma, 2.0));
}

double levlnorm(double d, double mu, double sigma, double order, int give_log)
{
  
  if (!R_FINITE(mu) ||
      !R_FINITE(sigma) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      sigma <= 0.0 ||
      d <= 0.0)
    return R_NaN;
  
  return exp(order * mu + 0.5 * R_pow(order, 2.0) * R_pow(sigma, 2.0)) * pnorm((log(d) - mu - order * R_pow(sigma, 2.0)) / sigma, 0, 1, 1, 0) + R_pow(d, order) * pnorm((log(d) - mu) / sigma, 0, 1, 0, 0);
}

