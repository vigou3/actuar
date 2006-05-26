/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the loglogistic distribution, and to simulate random
 *  variates. See ../R/llogis.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dllogis(double x, double shape, double scale, int give_log)
{
  double tmp1;
  double tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));

    tmp1 = shape * (log(x) - log(scale));
    tmp2 = R_pow(x / scale, shape);
    
    return  give_log ?
	log(shape) + tmp1 - log(x) - 2.0 * log(1.0 + exp(tmp1)) :
	shape * tmp2 / (x * R_pow(1.0 + tmp2, 2.0));
}

double pllogis(double x, double shape, double scale, int lower_tail, int log_p)
{

  double tmp1;
  double tmp2;

    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0) 
	error(_("invalid arguments"));

    if (x <= 0)
	return R_DT_0;

    tmp1 = shape + (log(x) - log(scale));
    tmp2 = R_pow(x / scale, shape);

    return (lower_tail ? R_D_exp(tmp1 - log(1.0 + exp(tmp1))):
	    R_D_exp(log(1.0 - exp(tmp1 - log(1.0 + exp(tmp1))))));
}

double qllogis(double p, double shape, double scale, int lower_tail, int log_p)
{

  double tmp1;
  double tmp2;

  if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 || 
	scale <= 0.0)
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, 1);

    tmp1 = 1.0 / shape;
    tmp2 = log(1-p) - log(p);

    return (lower_tail ? R_D_exp(log(scale) - tmp1 * tmp2) :
	    R_D_exp(log(scale) - tmp1 * (-tmp2)));
}


double rllogis(double shape, double scale)
{	
    double a;
	
    if (!R_FINITE(shape) ||
	!R_FINITE(scale) ||
	shape <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow(a / (1.0 - a), 1.0 / shape);
}
