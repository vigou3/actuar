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

double dtrbeta(double x, double shape1, double scale, double shape2, double shape3, int give_log)
{ 
    double tmp1, tmp2, tmp3;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	shape3 <= 0.0) 
	error(_("invalid arguments"));

    if (!R_FINITE(x) || x < 0.0) 
	return R_D_d0;

    tmp1 = x / scale;
    tmp2 = R_pow(tmp1, shape2);
    tmp3 = tmp2 / (1 + tmp2);
    
    /* !!! Modifier version log !!! */
    return  (give_log ?
	     dbeta(tmp3, shape3, shape1, 1) + log(shape2) - log(scale) + (shape2 - 1.0) * (log(x) - log(scale)) + 2.0 * (-log(1.0 + exp(shape2 * (log(x) - log(scale))))) :
	     dbeta(tmp3, shape3, shape1, 0) * shape2 * tmp2 / (scale * tmp1 * R_pow_di(1.0 + tmp2, 2)));
}

double ptrbeta(double q, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p)
{
    double u, tmp;

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	shape3 <= 0.0) 
	error(_("invalid arguments"));

    tmp = R_pow(q / scale, shape2);
    u = tmp / (1.0 + tmp);
    
    if (q <= 0)
	return R_DT_0;

    if (!R_FINITE(q))
        return 1;
    
    return (lower_tail ? R_D_exp(pbeta(u, shape3, shape1, 1, 1)):
	    R_D_exp(pbeta(u, shape3, shape1, 0, 1)));
}

double qtrbeta(double p, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p)
{
    double tmp;
    
    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	shape1 <= 0.0 || 
	scale <= 0.0 || 
	shape2 <= 0.0 || 
	shape3 <= 0.0) 
	error(_("invalid arguments"));
    
    R_Q_P01_boundaries(p, 0, R_PosInf);
    tmp = R_D_qIv(p);
    
    return  (lower_tail ? scale * R_pow(qbeta(tmp, shape3, shape1, 1, 0) / (1.0 - qbeta(tmp, shape3, shape1, 1 ,0)), 1.0 / shape2)  :
	     scale * R_pow((qbeta(tmp, shape3, shape1, 0, 0)) / (1.0 - qbeta(tmp, shape3, shape1, 0, 0)), 1.0 / shape2));
}

double rtrbeta(double shape1, double scale, double shape2, double shape3)
{
    double a;	

    if (!R_FINITE(shape1) ||
	!R_FINITE(scale) ||
	!R_FINITE(shape2) ||
	!R_FINITE(shape3) ||
	shape1 <= 0.0 ||
	scale <= 0.0 ||
	shape2 <= 0.0 ||
	shape3 <= 0.0)
	error(_("invalid arguments"));
	
    a = rbeta(shape3, shape1);

    return scale * R_pow(a / (1.0 - a), 1.0 / shape2);
}
