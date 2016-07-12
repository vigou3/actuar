/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Fonctions to compute density, cumulative distribution and quantile
 *  fonctions, raw and limited moments and to simulate random variates
 *  for the inverse Pareto distribution. See ../R/InversePareto.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "locale.h"
#include "dpq.h"

double dinvpareto(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u^shape * (1 - u) / x
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = x/scale.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape < 1) return R_PosInf;
	if (shape > 1) return ACT_D__0;
	/* else */
	return ACT_D_val(1.0 / scale);
    }

    tmp = log(x) - log(scale);
    logu = - log1pexp(-tmp);
    log1mu = - log1pexp(tmp);

    return ACT_D_exp(log(shape) + shape * logu + log1mu - log(x));
}

double pinvpareto(double q, double shape, double scale, int lower_tail,
                  int log_p)
{
    double u;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    if (q <= 0)
        return ACT_DT_0;

    u = exp(-log1pexp(log(scale) - log(q)));

    return ACT_DT_val(R_pow(u, shape));
}

double qinvpareto(double p, double shape, double scale, int lower_tail,
                  int log_p)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale / (R_pow(ACT_D_Lval(p), -1.0 / shape) - 1.0);
}

double rinvpareto(double shape, double scale)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    return scale / (R_pow(unif_rand(), -1.0 / shape) - 1.0);
}

double minvpareto(double order, double shape, double scale, int give_log)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape ||
        order >= 1.0)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(shape + order) * gammafn(1.0 - order)
        / gammafn(shape);
}

/* The function to integrate in the limited moment */
static void fn(double *x, int n, void *ex)
{
    int i;
    double *pars = (double *) ex, shape, scale, order;

    shape = pars[0];
    scale = pars[1];
    order = pars[2];

    for(i = 0; i < n; i++)
	x[i] = R_pow(x[i], shape + order - 1) * R_pow(1 - x[i], -order);
}

double levinvpareto(double limit, double shape, double scale, double order,
                    int give_log)
{
    double u;
    double ex[3], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    /* Parameters for the integral are pretty much fixed here */
    ex[0] = shape; ex[1] = scale; ex[2] = order;
    lower = 0.0; upper = limit / (limit + scale);
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;		     /* as instructed in WRE */
    iwork =   (int *) R_alloc(subdiv, sizeof(int));  /* idem */
    work = (double *) R_alloc(lenw, sizeof(double)); /* idem */

    Rdqags(fn, (void *) &ex,
	   &lower, &upper, &epsabs, &epsrel, &result,
	   &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    if (ier == 0)
    {
	u = exp(-log1pexp(log(scale) - log(limit)));
	return R_pow(scale, order) * shape * result
	    + ACT_DLIM__0(limit, order) * (0.5 - R_pow(u, shape) + 0.5);
    }
    else
	error(_("integration failed"));
}
