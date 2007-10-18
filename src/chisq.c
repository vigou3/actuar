/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the Chi-square
 *  distribution. See ../R/ChisqSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mchisq(double order, double df, double ncp, int give_log)
{
    if (!R_FINITE(df) ||
	!R_FINITE(ncp) ||
	!R_FINITE(order) ||
	df <= 0.0 ||
	ncp < 0.0 ||
	order <= -df/2)
	return R_NaN;

	//centred chi-square distribution
	if(ncp == 0) 
		return R_pow(2, order) * gammafn(order + df/2) / gammafn(df/2);
		
	//non centred chi-square distribution
	if(order >= 1.0)
	{
		int i,j = 0; //loop index	
		int	n = order;

		double Fact[n+1];//array with 0!, 1!, 2!, ... n!	
		/* init */ 
		Fact[0] = 1; //0!      
		for( i=1; i< n+1 ; i++) Fact[i]= i * Fact[i-1];
	
		double Moment[n+1];	//array with 1, E(X), E(X^2), E(X^3),... E(X^n)
							//i.e. E(X^i) = Moment[i]
		/* init */
		Moment[0] = 1;
		Moment[1] = df + ncp; //E(X)
		for( i=2; i< n+1 ; i++)
		{ 
			Moment[i] = R_pow(2, i-1)*Fact[i-1]*(df + i*ncp);			
			for( j=1; j<i ; j++)
				Moment[i] += Fact[i-1]/Fact[i-j]*R_pow(2, j-1)*(df + j*ncp)*Moment[i-j];
		}	
			
		return Moment[n];	
	}			
}

double levchisq(double limit, double df, double ncp, double order, int give_log)
{
    double u, tmp;

    if (!R_FINITE(df) ||
	!R_FINITE(ncp) ||
	!R_FINITE(order) ||
	df <= 0.0 ||
	ncp < 0.0 ||
	order <= -df/2)
	return R_NaN;

    if (limit <= 0.0)
	return 0;

	if(ncp == 0)
	{
		tmp = order + df/2;

		u = exp(log(limit) - log(2));

		return R_pow(2, order) * gammafn(tmp) *
		pgamma(u, tmp, 1.0, 1, 0) / gammafn(df/2) +
		R_VG__0(limit, order) * pgamma(u, df/2, 1.0, 0, 0);
	}
	else
		return 0.0;
}

double mgfchisq(double x, double df, double ncp, int give_log)
{
	/*check arguments */
	if (!R_FINITE(df) ||
	    !R_FINITE(ncp) ||
	    df <= 0.0 ||
	    ncp < 0.0 ||
	    2 * x > 1.)
	  return R_NaN;
	  
	
	if(x == 0.0)
	  return R_D_exp(0.0);	
	
	return	R_D_exp( ncp*x / (1. - 2*x) - df/2 * log(1. - 2*x) ) ;
}
