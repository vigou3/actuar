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

double mnorm(int order, double mean, double sd, int give_log)
{
	if (!R_FINITE(mean) || 
		!R_FINITE(sd) ||	
		!R_FINITE(order) ||
		sd <= 0.0)
		return R_NaN;
		
	if(order == 0.0)
		return 1.0;
		
	int i = 0; //loop index	
	int length = order;

	double Fact[length+1];//array with 0!, 1!, 2!, ... n!
	/* init */ 
	Fact[0] = 1;       
	for( i=1; i< length+1 ; i++) Fact[i]= i * Fact[i-1];		
		
	double res = 0;
	for( i=0; i< length/2+1; i++) res += Fact[length] / (pow(2,i) * Fact[i] * Fact[length-2*i] ) * pow(sd,2*i) * pow(mean,order-2*i);
						
	return res;													
}

double mgfnorm(double t, double mean, double sd, int give_log)
{
	/*check arguments */
	if (!R_FINITE(mean) ||
	    !R_FINITE(sd) ||
	    sd <= 0.0 )
	  return R_NaN;
	  
	
	if(t == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( t*mean + 0.5*t*t*sd*sd ) ;
}
