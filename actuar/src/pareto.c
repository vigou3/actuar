#include <R.h>
#include <Rmath.h>

void dpareto(double *x, double *alpha, double *lambda, double *y, int *length)
{
    int i;

    if (R_IsNaN(*x) || R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("x, alpha et lambda doivent être des nombres");
    if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
    if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");

/*     if (*x < 0) return *y = 0; */

    for(i = 0; i <= *length; i++)
	y[i] = exp(log(*alpha) + *alpha*log(*lambda) - (*alpha + 1.0)*log(x[i] + *lambda));
}


void ppareto(double *x, double *alpha, double *lambda, double *y, int *length)
{
    int i;
    
    if (R_IsNaN(*x) || R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("x, alpha et lambda doivent être des nombres");
    if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
    if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");
    
/*     if (*x < 0) return *y = 0; */
    
    for (i = 0; i <= *length; i++)
	y[i] = 1.0 - exp(*alpha * (log(*lambda) - log(x[i] + *lambda)));
}


void qpareto(double *p, double *alpha, double *lambda, double *y, int *length)
{
    int i;
    
    if (R_IsNaN(*p) || R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("x, alpha et lambda doivent être des nombres");
    if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
    if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");
    if (*p < 0) Rf_error ("x doit être compris entre 0 et 1");
    if (*p > 1) Rf_error ("x doit être compris entre 0 et 1");
    
    for (i = 0; i <= *length; i++)
	y[i] = exp(log(*lambda) + log((1.0 - pow(1.0 - p[i], 1.0 / *alpha))/pow(1.0 - p[i], 1.0 / *alpha)));
}


void rpareto(int *n, double *alpha, double *lambda, double *y)
{
    int i;

    if (R_IsNaN(*alpha) || R_IsNaN(*lambda)) Rf_error ("alpha et lambda doivent être des nombres");
    if (*alpha <= 0)  Rf_error ("alpha et lambda doivent être positifs");
    if (*lambda <= 0) Rf_error ("alpha et lambda doivent être positifs");
    
    GetRNGstate();
    
    for (i = 0; i <= *n; i++)
	y[i] = *lambda * pow(unif_rand(), -1.0 / *alpha) - *lambda;

    PutRNGstate();
}
