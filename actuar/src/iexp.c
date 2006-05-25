/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions of the Pareto distribution, and to simulate random
 *  variates. See ../R/pareto.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double diexp(double x, double scale, int give_log)
{
    if (!R_FINITE(scale) ||
	scale <= 0.0 || 
	x < 0.0) 
	error(_("invalid arguments"));
    
    return  give_log ?
	log(scale) - scale / x - 2.0 * log (x) :
	scale * exp(-scale / x) / R_pow(x, 2.0);
    
}

double riexp(double scale)
{
    if (!R_FINITE(scale) ||
	scale <= 0.0)
	error(_("invalid arguments"));

    return scale / log(1.0 / unif_rand());
}
