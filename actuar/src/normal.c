/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Gamma
 *  distribution. See ../R/NormalSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"


double mgfnorm(double x, double mean, double sd, int give_log)
{
	/*check arguments */
	if (!R_FINITE(mean) ||
	    !R_FINITE(sd) ||
	    sd <= 0.0 )
	  return R_NaN;
	  
	
	if(x == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( x*mean + 0.5*x*x*sd*sd ) ;
}
