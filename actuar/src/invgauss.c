/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse gamma distribution, and to simulate random
 *  variates. See ../R/invgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgauss(double x, double mean, double scale, int give_log)
{
	double tmp;

    if (!R_FINITE(mean) ||
	!R_FINITE(scale)  ||
	mean <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

     if (!R_FINITE(x) || x < 0.0) 
      return R_D_d0;

	tmp = (x - mean) / mean;
    
    return  give_log ?
      (0.5) * (log(scale) - log(2.0 * M_PI * x * x * x)) - scale * tmp * tmp / (2.0 * x) :
      R_pow(scale, 0.5) * exp(-scale * tmp * tmp / (2.0 * x)) / R_pow(2.0 * M_PI * x * x * x, 0.5);
}

double pinvgauss(double q, double mean, double scale, int lower_tail, int log_p)
{
  double tmp, tmp1, tmp2;
  
    if (!R_FINITE(mean) || 
	!R_FINITE(scale) ||
	mean <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = (q - mean) / mean;
    tmp1 = (q + mean) / mean;
    tmp2 = scale / q;
  
      return (lower_tail ?
        pnorm(tmp * R_pow(tmp2, 0.5), 0, 1, 1, 0) + exp(2.0 * scale / mean) * pnorm(-tmp1 * R_pow(tmp2, 0.5), 0, 1, 1, 0):
	      1.0 - (pnorm(tmp * R_pow(tmp2, 0.5), 0, 1, 1, 0) + exp(2.0 * scale / mean) * pnorm(-tmp1 * R_pow(tmp2, 0.5), 0, 1, 1, 0)));
}

double qinvgauss(double p, double mean, double scale, int lower_tail, int log_p)
{
  double tmp, Y, r1, test;

  if (!R_FINITE(mean) || 
	!R_FINITE(scale) ||
	mean <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);

  Y = qchisq(tmp,1, 1, 0); 
  r1 = mean * (2.0 * scale + mean * Y - R_pow(4.0 * scale * mean * Y + mean * mean * Y * Y, 0.5) / (2.0 * scale));
  test = mean / (mean + r1);

    return (lower_tail ? r1 :
	    1.0 - r1);
}

double rinvgauss(double mean, double scale)
{

  double a, r1, test;

    if (!R_FINITE(mean) ||
	!R_FINITE(scale) ||
	mean <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = rchisq(1);
    r1 = mean + (mean / (2.0 * scale)) * (mean * a - R_pow(mean * a * (4.0 * scale + mean * a), 0.5));
    test = rbinom(1, mean / (mean + r1));

    return  test * r1 + (1.0 - test) * mean * (mean / r1);
}
