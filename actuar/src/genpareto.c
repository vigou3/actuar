/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the generalized Pareto distribution, and to simulate random
 *  variates. See ../R/genpareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dgenpareto(double x, double shape1, double scale, double shape2, int give_log)
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

    tmp = x / (x + scale);
    
    return  give_log ?
      dbeta(tmp, shape2, shape1, 1) + log(scale) - 2.0 * log(x + scale) :
      dbeta(tmp, shape2, shape1, 0) * scale / R_pow(x + scale, 2.0);
}

double pgenpareto(double q, double shape1, double scale, double shape2, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    u = q / (q + scale);

    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
        return 1;
    
    return (lower_tail ? R_D_exp(pbeta(u, shape2, shape1, 1, 1)):
	    R_D_exp(pbeta(u, shape2, shape1, 0, 1)));
}

double qgenpareto(double p, double shape1, double scale, double shape2, int lower_tail, int log_p)
{

  double tmp;

  if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0) 
	error(_("invalid arguments"));

    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);

    return (lower_tail ? scale * (qbeta(tmp, shape2, shape1, 1, 0) / (1.0 - qbeta(tmp, shape2, shape1, 1, 0))) :
	    scale * qbeta(tmp, shape2, shape1, 0, 0) / (1.0 - qbeta(tmp, shape2, shape1, 0, 0)));
}

double rgenpareto(double shape1, double scale, double shape2)
{	
    double a;
	
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0)
	error(_("invalid arguments"));

    a = rbeta(shape2, shape1);

    return scale * a / (1.0 - a);
}
