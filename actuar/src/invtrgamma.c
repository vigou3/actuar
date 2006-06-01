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

double dinvtrgamma(double x, double shape1, double scale, double shape2, int give_log)
{
	double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale)  ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	x < 0.0) 
	error(_("invalid arguments"));

	tmp = R_pow(1.0 / x, shape2);
    
    return  give_log ?
      log(shape2) - (shape2 + 1.0) * log(x) + dgamma(tmp, shape1, 1.0 / scale, 1) :
      shape2 *  R_pow(x, -shape2 - 1.0) * dgamma(tmp, shape1, 1.0 / scale, 0);
}

double pinvtrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 ||
        shape2 <= 0.0) 
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = R_pow(1.0 / q, shape2);
    
    return (lower_tail ? R_D_exp(pgamma(tmp, shape1, 1.0 / scale, 0,1)):
	    R_D_exp(pgamma(tmp, shape1, 1.0 / scale, 1,1)));
}

double qinvtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 ||
        shape2 <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, 1);
  tmp = R_D_qIv(p);

    return (lower_tail ? R_pow(qgamma(tmp, shape1, 1.0 / scale, 0, 0), -1.0 / shape2) :
	    R_pow(qgamma(tmp, shape1, 1.0 / scale, 1, 0), -1.0 / shape2));
}

double rinvtrgamma(double shape1, double scale, double shape2)
{
  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 ||
      shape2 <= 0.0)
	error(_("invalid arguments"));

    return qgamma(unif_rand(), shape1, 1.0 / scale, 1, 0);
}
