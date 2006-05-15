#include "R.h"
#include "Rmath.h"

void qpareto(double *x, double *alpha, double *lambda, double *y, int *length)
{
	int i;

	if (R_IsNaN(*x) || R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("x, alpha et lambda doivent être des nombres");
	if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
	if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");
	if (*x < 0) Rf_error ("x doit être compris entre 0 et 1");
	if (*x > 1) Rf_error ("x doit être compris entre 0 et 1");
	
	for (i = 0; i <= *length; i++)
		{
			y[i] = exp(log(*lambda) + log((1 - pow(1-x[i], 1 / *alpha))/pow(1-x[i],1 / *alpha)));
		}
}
