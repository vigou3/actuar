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
		return (R_pow(max,order+1) - R_pow(min,order+1)) /( (max - min)*(order+1) );		
}

double levunif(double limit, double min, double max, double order, int give_log)
{
	/*check arguments */
	if (!R_FINITE(min) || !R_FINITE(max) || min >= max)
		return R_NaN;
	
	double tmp, res;
	
	if(limit <= min)
		return R_pow(limit,order);
		
	if(limit >= max)
		return munif(order, min, max, give_log);
	
	if(order == -1.0)
	{
		tmp = (log(abs(limit)) - log(abs(min))) / (max - min);		
		res = tmp + (max - limit) / (limit*(max - min));
	}
	else
	{
		tmp = (R_pow(limit,order+1) - R_pow(min,order+1)) /( (max - min)*(order+1) );
		res = tmp + R_pow(limit, order)*(max - limit)/(max - min);
	}
	return res;			
}

double mgfunif(double x, double min, double max, int give_log)
{	
	/*check arguments */
	if (!R_FINITE(min) || !R_FINITE(max) || min >= max)
		return R_NaN;	  
	
	if(x == 0.0)
		return R_D_exp(0.0);	
	
	double tmp1 = exp(x*max)-exp(x*min);
	double tmp2 = x*(max-min);
	
	//we can't use the macro R_D_exp since the log
	//of the mgf never exists for t<0
	if(give_log)
		return	log(tmp1) - log(tmp2) ;
	else
		return tmp1/tmp2 ;
}
