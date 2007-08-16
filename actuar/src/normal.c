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

double mnorm(double order, double mean, double sd, int give_log)
{
	printf("param %f,%f,%f\n",order,mean,sd);
	
	if (!R_FINITE(mean) || 
		!R_FINITE(sd) ||	
		!R_FINITE(order) ||
		sd <= 0.0)
		return R_NaN;
		
	if(order == 0.0)
		return 1.0;
		
	int i = 0; //loop index	
	int n = order;

	double Fact[n+1];//array with 0!, 1!, 2!, ... n!
	/* init */ 
	Fact[0] = 1;       
	for( i=1; i< n+1 ; i++) Fact[i]= i * Fact[i-1];		
		
	double res = 0;
	for( i=0; i< n/2+1; i++) 
		res += Fact[n] / (pow(2,i) * Fact[i] * Fact[n-2*i] ) * pow(sd,2*i) * pow(mean,n-2*i);
						
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
