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
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    if (!R_FINITE(x)  ||
	x < 0.0) 
      return R_D_d0;

    tmp = R_pow(x / scale, shape2);
    
    return  (give_log ?
	log(shape1) + log(shape2) + shape2 * (log(x) - log(scale)) - log(x) - (shape1 + 1.0) * (log(R_pow(scale, shape2) + R_pow(x, shape2)) - shape2 * log(scale)) :
	shape1 * shape2 * tmp / (x * R_pow(1.0 + tmp, shape1 + 1.0)));
}

double pburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  double tmp1, tmp2;

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

  double tmp, tmp1, tmp2;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
        shape2 <= 0.0)
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, R_PosInf );
    tmp =  R_D_qIv(p);

    tmp1 = 1.0 / shape1;
    tmp2 = 1.0 / shape2;

    return (lower_tail ? scale * R_pow((1.0 - R_pow(1.0 - tmp, tmp1)) / R_pow(1.0 - tmp, tmp1), tmp2) :
	    scale * R_pow((1.0 - R_pow(tmp, tmp1)) / R_pow(tmp, tmp1), tmp2));
}

double rburr(double shape1, double scale, double shape2)
{
    double a, tmp;
	
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

double mburr(double k, double shape1, double scale, double shape2, int give_log)
{
	
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(k) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	k <= -shape2 ||
	k >= shape1 * shape2)
	error(_("invalid arguments"));
	
    return R_pow(scale, k) * gammafn(1.0 + k / shape2) * gammafn(shape1 - k / shape2) / gammafn(shape1);
}

double levburr(double x, double shape1, double scale, double shape2, double order, int give_log)
{
  double u;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(x) ||
	!R_FINITE(order) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	x <= 0.0 ||
	order <= -shape2)
	error(_("invalid arguments"));

    u = 1.0 / (1.0 + R_pow(x / scale, shape2));

    return R_pow(scale, order) * gammafn(1.0 + order / shape2) * gammafn(shape1 - order / shape2) * pbeta(1.0 - u, 1.0 + order / shape2, shape1 - order / shape2, 1, 0) / gammafn(shape1) + R_pow(x, order) * R_pow(u, shape1);
}
