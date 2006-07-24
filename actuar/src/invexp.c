/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse exponential distribution, to calculate raw moments and limited moments 
 *  of the random variable and to simulate random variates. See ../R/invexp.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvexp(double x, double scale, int give_log)
{
  
  double tmp;
  
  if (!R_FINITE(scale) || scale <= 0.0) 
    return R_NaN;
  
  if (!R_FINITE(x) || x < 0.0) 
    return R_D_d0;
  
  tmp = scale / x;
  
  return  give_log ?
    log(scale) - tmp - 2.0 * log (x) :
    scale * exp(-tmp) / R_pow(x, 2.0);
}

double pinvexp(double q, double scale, int lower_tail, int log_p)
{
  double tmp;
  
  if (!R_FINITE(scale) || scale <= 0.0)
    return R_NaN;
  
  if (q <= 0)
    return R_DT_0;
  
  tmp = scale / q;
  
  return (lower_tail ? R_D_exp(-tmp):
	  R_D_exp(log(1.0 - exp(-tmp))));
}

double qinvexp(double p, double scale, int lower_tail, int log_p)
{
  double tmp;
  
  if (!R_FINITE(scale) || scale <= 0.0)
    return R_NaN;
  
  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);
  
  return (lower_tail ? -scale / log(tmp) :
	  -scale / log(1.0 - tmp)) ;
}


double rinvexp(double scale)
{
  if (!R_FINITE(scale) || scale <= 0.0)
    return R_NaN;
  
  return scale / log(1.0 / unif_rand());
}

double minvexp(double k, double scale, int give_log)
{
  if (!R_FINITE(scale) || !R_FINITE(k) || scale <= 0.0 || k >= 1.0)
    return R_NaN;
  
  return R_pow(scale, k) * gammafn(1.0 - k);
}

double levinvexp(double d, double scale, double order, int give_log)
{
  if (!R_FINITE(scale) || !R_FINITE(d) || R_FINITE(order) || scale <= 0.0 || order <= 0.0 || order >= 1.0)
    return R_NaN;
  
  return R_pow(scale, order) * gammafn(1.0 - order) * (1.0 - pgamma(1.0 / d, 1.0 - order, 1.0 / scale, 1, 0)) + R_pow(d, order) * (1.0 - exp(-scale / d));
}
