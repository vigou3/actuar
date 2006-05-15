#include "R.h"
#include "Rmath.h"
#include "dpq.h"

void rpareto(double *alpha, double *lambda, double *y, int *length)
{
	int i;

	if (R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("alpha et lambda doivent être des nombres");
	if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
	if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");

	for (i=0; i<= *length; i++)
		{
			GetRNGstate();
			y[i] = *lambda * pow(1.0-unif_rand(),-1.0 / *alpha) - *lambda;
			PutRNGstate();
		
		}

	

}
