/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the Pareto distribution, and to simulate random
 *  variates. See ../R/lgamma.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dlgamma(double x, double shapelog, double ratelog, int give_log)
{
    if (!R_FINITE(shapelog) || 
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 || 
	ratelog <= 0.0) 
	error(_("invalid arguments"));

      if (!R_FINITE(x) || x < 0.0) 
      return R_D_d0;
    
    return (give_log ?
	-log(x) + dgamma(log(x), shapelog, 1.0 / ratelog, 1)  :
	dgamma(log(x), shapelog, 1.0 / ratelog, 0) / x);
}

double plgamma(double q, double shapelog, double ratelog, int lower_tail, int log_p)
{
  
    if (!R_FINITE(shapelog) || 
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 || 
	ratelog <= 0.0)
	error(_("invalid arguments"));

    if (q <= 0)
	return R_DT_0;
    
    return (lower_tail ? R_D_exp(pgamma(log(q), shapelog, 1.0 / ratelog, 1, 1)):
	    R_D_exp(pgamma(log(q), shapelog, 1.0 / ratelog, 0, 1)));
}

double qlgamma(double p, double shapelog, double ratelog, int lower_tail, int log_p)
{
  double tmp;

  if (!R_FINITE(shapelog) || 
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 || 
	ratelog <= 0.0)
	error(_("invalid arguments"));

  R_Q_P01_boundaries(p, 1, R_PosInf);
  tmp = R_D_qIv(p);

    return (lower_tail ? exp(qgamma(tmp, shapelog, 1.0 / ratelog, 1, 0)) :
	    exp(qgamma(tmp, shapelog, 1.0 / ratelog, 0, 0)));
}

double rlgamma(double shapelog, double ratelog)
{
  double a;

    if (!R_FINITE(shapelog) || 
	!R_FINITE(ratelog) ||
	shapelog <= 0.0 || 
	ratelog <= 0.0)
	error(_("invalid arguments"));

    a = rgamma(shapelog, 1.0 / ratelog);

    return exp(a);
}
