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

	//TODO
	
	return 0.0;		
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

double mgfinvGauss(double t, double nu, double lambda, int give_log)
{
	/*check arguments */
	if (!R_FINITE(nu) ||
	    !R_FINITE(lambda) ||
	    nu <= 0.0 ||
	    lambda < 0.0 ||
		t > lambda/(2*nu*nu) )
	  return R_NaN;
	  
	
	if(t == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( lambda / nu *( 1 - sqrt( 1-2*nu*nu*t/lambda ) ) ) ;
}
