/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the inverse exponential distribution, and to simulate random
 *  variates. See ../R/invexp.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvexp(double x, double scale, int give_log)
{

  double tmp;

    if (!R_FINITE(scale) ||
	scale <= 0.0) 
	error(_("invalid arguments"));

     if (!R_FINITE(x)  ||
	x < 0.0) 
      return R_D_d0;

    tmp = scale / x;
    
    return  give_log ?
	log(scale) - tmp - 2.0 * log (x) :
	scale * exp(-tmp) / R_pow(x, 2.0);
    
}

double pinvexp(double q, double scale, int lower_tail, int log_p)
{
  double tmp;
  
    if (!R_FINITE(scale) ||
	scale <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;

    tmp = scale / q;
    
    return (lower_tail ? R_D_exp(-tmp):
	    R_D_exp(log(1.0 - exp(-tmp))));
}

double qinvexp(double p, double scale, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(scale) || 
	scale <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 0, R_PosInf);
  tmp = R_D_qIv(p);

    return (lower_tail ? -scale / log(tmp) :
	   -scale / log(1.0 - tmp)) ;
}


double rinvexp(double scale)
{
    if (!R_FINITE(scale) ||
	scale <= 0.0)
	error(_("invalid arguments"));

    return scale / log(1.0 / unif_rand());
}
