/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Gamma
 *  distribution. See ../R/UniformSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"


double mgfunif(double x, double min, double max, int give_log)
{
	/*check arguments */
	if (!R_FINITE(min) || !R_FINITE(max) || min > max)
		return R_NaN;
	  
	
	if(x == 0.0)
		return R_D_exp(0.0);	
	
	double tmp1 = exp(x*max)-exp(x*min);
	double tmp2 = x*(max-min);
	
	return	R_D_exp( log(tmp1) - log(tmp2) ) ;
}
