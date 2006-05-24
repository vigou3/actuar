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

double rinverseparalogistic(double tau, double scale)
{
    double a;
	
    if (!R_FINITE(tau) ||
	!R_FINITE(scale) ||
	tau <= 0.0 ||
	scale <= 0.0)
	error(_("invalid arguments"));

    a = unif_rand();

    return scale * R_pow((R_pow(a, 1.0 / tau)) / (1.0 - R_pow(a, 1.0 / tau)), 1.0 / tau);
}
