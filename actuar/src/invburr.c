/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse Burr distribution, and to simulate random
 *  variates. See ../R/invburr.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvburr(double x, double shape1, double scale, double shape2, int give_log)
{
  double tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    tmp = log(x) - log(scale);
    
    return  give_log ?
	log(shape1) + log(shape2) + shape2 * shape1 * tmp - log(x) - (shape1 + 1.0) * log(1.0 + exp(shape2 * tmp)) :
	shape1 * shape2 * R_pow(x / scale, shape2 * shape1) / (x * R_pow(1 + R_pow(x / scale, shape2), shape1 + 1.0));
    
}

double pinvburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
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

    tmp = log(q) - log(scale);

    return (lower_tail ? R_D_exp(shape1 * shape2 * tmp - shape1 * log(1.0 + exp(shape2 * tmp))):
	    R_D_exp(log(1.0 - exp(shape1 * shape2 * tmp - shape1 * log(1.0 + exp(shape2 * tmp))))));
}

double qinvburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
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

    R_Q_P01_boundaries(p, 0, 1);
    tmp = R_D_qIv(p);

    tmp1 = R_pow(tmp, 1.0 / shape1);

    return (lower_tail ? scale * R_pow(tmp1 / (1.0 - tmp1), 1.0 / shape2) :
	    scale * R_pow((1.0 - tmp1) / tmp1, 1.0 / shape2));
}


double rinvburr(double shape1, double scale, double shape2)
{	
    double a;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow(R_pow(a, 1.0 / shape1) / (1.0 - R_pow(a, 1.0 / shape1)), 1.0 / shape2);
}
