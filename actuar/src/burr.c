/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the Burr distribution, and to simulate random
 *  variates. See ../R/burr.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dburr(double x, double shape1, double scale, double shape2, int give_log)
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

/*  tmp = x / scale;
 */
   
    tmp = R_pow(x / scale, shape2);
    
/*   return  (give_log ? 
 * 	     dbeta(x, 1.0, shape1, 1) + log(shape2) - log(scale) + (shape2 - 1.0)*(log(x) - log(scale)) + 2.0 * (-log(1.0 + exp(shape2 * (log(x) - log(scale))) *)) : 
 * 	     dbeta(x, 1.0, shape1, 0) * (shape2 / scale) * R_pow(tmp, shape2 - 1.0) * R_pow(1.0 / (1.0 + R_pow(tmp), shape2), 2.0)); 
 */

    return  (give_log ?
	log(shape1) + log(shape2) + shape2 * (log(x) - log(scale)) - log(x) - (shape1 + 1.0) * (log(R_pow(scale, shape2) + R_pow(x, shape2)) - shape2 * log(scale)) :
	shape1 * shape2 * tmp / (x * R_pow(1.0 + tmp, shape1 + 1.0)));

}

double pburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  double tmp1;
  double tmp2;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp1 = R_pow(scale, shape2);
    tmp2 = R_pow(q, shape2);

    return (lower_tail ? R_D_exp(log(1.0 - exp(shape1 * (shape2 * log(scale) - log(tmp1 + tmp2))))):
	    R_D_exp(shape1 * (shape2 * log(scale) - log(tmp1 + tmp2))));
}

double qburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  double tmp;
  double tmp1;
  double tmp2;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
      shape2 <= 0.0)
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);
    tmp =  R_D_qIv(p);

    tmp1 = 1.0 / shape1;
    tmp2 = 1.0 / shape2;

    return (lower_tail ? scale * R_pow((1.0 - R_pow(1.0 - tmp, tmp1)) / R_pow(1.0 - tmp, tmp1), tmp2) :
	    scale * R_pow((1.0 - R_pow(tmp, tmp1)) / R_pow(tmp, tmp1), tmp2));
}

double rburr(double shape1, double scale, double shape2)
{
    double a;
    double tmp;
	
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();
    tmp = 1.0 / shape1;
	
    return scale * R_pow((1.0 - R_pow(a, tmp)) / (R_pow(a, tmp) ), 1.0 / shape2);
}
