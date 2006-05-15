#include "R.h"
#include "Rmath.h"

void ppareto(double *x, double *alpha, double *lambda, double *y, int *length)
{

	int i;

	if (R_IsNaN(*x) || R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("x, alpha et lambda doivent être des nombres");
	if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
	if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");

	if (*x < 0) return *y = 0;

	for (i=0; i <= *length; i++)
		{
			y[i] = 1 - exp(*alpha*(log(*lambda) - log(x[i] + *lambda)));
		}
}
