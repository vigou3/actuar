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

double dinvgamma(double x, double shape, double scale, int give_log)
{
  double tmp1, tmp2;

  if (!R_FINITE(shape) ||
      !R_FINITE(scale)  ||
      shape <= 0.0 ||
      scale <= 0.0)
    return R_NaN;

  if (!R_FINITE(x) || x < 0.0)
    return R_D_d0;

  tmp1 = 1.0 / x;
  tmp2 = 1.0 / scale;

  return  give_log ?
    -2.0 * log(x) + dgamma(tmp1, shape, tmp2, 1) :
    R_pow(x, -2.0) * dgamma(tmp1, shape, tmp2, 0);
}

double pinvgamma(double q, double shape, double scale, int lower_tail, int log_p)
{
  double tmp1, tmp2;

  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      shape <= 0.0 ||
      scale <= 0.0)
    return R_NaN;;

  if (q <= 0)
    return R_DT_0;

  tmp1 = 1.0 / q;
  tmp2 = 1.0 / scale;

  return (lower_tail ? R_D_exp(pgamma(tmp1, shape, tmp2, 0,1)):
	  R_D_exp(pgamma(tmp1, shape, tmp2, 1,1)));
}

double qinvgamma(double p, double shape, double scale, int lower_tail, int log_p)
{
  double tmp1, tmp2;

  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      shape <= 0.0 ||
      scale <= 0.0)
    return R_NaN;;

  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp1 = R_D_qIv(p);
  tmp2 = 1.0 / scale;

  return (lower_tail ? 1.0 / qgamma(tmp1, shape, tmp2, 0, 0) :
	  1.0 / qgamma(tmp1, shape, tmp2, 1, 1));
}

double rinvgamma(double shape, double scale)
{
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      shape <= 0.0 ||
      scale <= 0.0)
    return R_NaN;;

  /*    return 1.0 / qgamma(unif_rand(), shape, 1.0 / scale, 1, 0); */
  return 1.0 / rgamma(shape, 1.0 / scale);
}

double minvgamma(double k, double shape, double scale, int give_log)
{
  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      !R_FINITE(k) ||
      shape <= 0.0 ||
      scale <= 0.0 ||
      k >= shape)
    return R_NaN;;

  return R_pow(scale, k) * gammafn(shape - k) / gammafn(shape);
}

double levinvgamma(double d, double shape, double scale, double order, int give_log)
{
  double tmp2;

  if (!R_FINITE(shape) ||
      !R_FINITE(scale) ||
      !R_FINITE(d) ||
      !R_FINITE(order) ||
      shape <= 0.0 ||
      scale <= 0.0 ||
      order >= shape ||
      d <= 0.0)
    return R_NaN;;
  tmp2 = 1.0 / scale;

  return R_pow(scale, order) * gammafn(shape - order) * (pgamma(1.0 / d, shape - order, tmp2, 0, 0)) / gammafn(shape) + R_pow(d, order) * pgamma(1.0 / d, shape, tmp2, 1, 0);
}
