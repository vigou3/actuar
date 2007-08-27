/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the inverse gaussian
 *  distribution. See ../R/InvGaussSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double minvGauss(double order, double nu, double lambda, int give_log)
{
    if (!R_FINITE(nu) ||
	!R_FINITE(lambda) ||
	!R_FINITE(order) ||
	nu <= 0.0 ||
	lambda <= 0.0 )
	return R_NaN;

	int n = order;
	int i = 0; //loop index
	
	if(n == 0)
		return 0.0;		
	
	double Fact[2*n-1];//array with 0!, 1!, 2!, ... n!	
	/* init */ 
	Fact[0] = 1; //0!      
	for( i=1; i< 2*n+1 ; i++) Fact[i]= i * Fact[i-1];
	
	double res = 0;
	
	for(i=0; i < n; i++) 
		res += R_pow(nu,n) * Fact[n-1+i] / (Fact[i]*Fact[n-1-i]) * R_pow(2*lambda/nu, -i);
						
	return res;		
}

double levinvGauss(double limit, double nu, double lambda, double order, int give_log)
{
    double z, y, tmpZ, tmpY;

    if (!R_FINITE(nu) ||
	!R_FINITE(lambda) ||
	!R_FINITE(order) ||
	nu <= 0.0 ||
	lambda < 0.0 )
	return R_NaN;

    if (limit <= 0.0)
	return 0;

	// when order = 1
	z = (limit - nu)/nu;
	y = (limit + nu)/nu;
	tmpZ = pnorm( z * sqrt(lambda / limit), 0., 1.0, 0, 0);
	tmpY = pnorm( -y * sqrt(lambda / limit), 0., 1.0, 0, 0);
	
	return limit - nu * z * tmpZ - nu * y * exp( 2*lambda/nu ) * tmpY;
}

double mgfinvGauss(double x, double nu, double lambda, int give_log)
{
	/*check arguments */
	if (!R_FINITE(nu) ||
	    !R_FINITE(lambda) ||
	    nu <= 0.0 ||
	    lambda < 0.0 ||
		x > lambda/(2*nu*nu) )
	  return R_NaN;
	  
	
	if(x == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( lambda / nu *( 1 - sqrt( 1-2*nu*nu*x/lambda ) ) ) ;
}
