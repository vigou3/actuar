/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the transformed beta distribution, and to simulate random
 *  variates. See ../R/trbeta.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dtrgamma(double x, double shape1, double scale, double shape2, int give_log)
{ 
  double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    if (!R_FINITE(x)  ||
	x < 0.0) 
      return R_D_d0;

    tmp = R_pow(x, shape2);
    
    return  (give_log ? dgamma(tmp, shape1, 1.0 / scale, 1) + log(shape2) + (shape2 - 1.0) * log(x) :
	     shape2 * R_pow(x, shape2 - 1.0) * dgamma(tmp, shape1, 1.0 / scale, 0));
}

double ptrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
    double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0)
	error(_("invalid arguments"));

    tmp = R_pow(q, shape2);

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
        return 1;
    
    return (lower_tail ? R_D_exp(pgamma(tmp, shape1, 1.0 / scale, 1, 1)):
	    R_D_exp(pgamma(tmp, shape1, 1.0 / scale, 0, 1)));
}

double qtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
  double tmp;
  double tmp1;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);
    tmp1 = 1.0 / shape2;

    return (lower_tail ? R_pow(qgamma(tmp, shape1, 1.0 / scale, 1, 0), tmp1) :
	    R_pow(qgamma(tmp, shape1, 1.0 / scale, 0, 0), tmp1));
}

double rtrgamma(double shape1, double scale, double shape2)
{
    double a;	

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	error(_("invalid arguments"));
	
    a = rgamma(shape1, 1.0 / scale);

    return R_pow(a, 1.0 / shape2);
}
