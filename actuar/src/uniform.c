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

double munif(double order, double min, double max, int give_log)
{
	/*check arguments */
	if (!R_FINITE(min) || !R_FINITE(max) || min >= max)
		return R_NaN;
	
	if(order == -1.0)
		return (log(abs(max)) - log(abs(min))) / (max - min);
	else
		return (pow(max,order+1) - pow(min,order+1)) /( (max - min)*(order+1) );		
}


double mgfunif(double t, double min, double max, int give_log)
{	
	/*check arguments */
	if (!R_FINITE(min) || !R_FINITE(max) || min >= max)
		return R_NaN;	  
	
	if(t == 0.0)
		return R_D_exp(0.0);	
	
	double tmp1 = exp(t*max)-exp(t*min);
	double tmp2 = t*(max-min);
	
	//we can't use the macro R_D_exp since the log
	//of the mgf never exists for t<0
	if(give_log)
		return	log(tmp1) - log(tmp2) ;
	else
		return tmp1/tmp2 ;
}
