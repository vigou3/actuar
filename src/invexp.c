/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the inverse exponential distribution. See ../R/InverseExponential.R
 *  for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dinvexp(double x, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  u * e^(-u) / x
     *
     *  with u = scale/x.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(scale))
	return x + scale;
#endif
    if (!R_FINITE(scale) || scale <= 0.0)
        return R_NaN;

    /* handle also x == 0 here */
    if (!R_FINITE(x) || x <= 0.0)
        return ACT_D__0;

    double logu = log(scale) - log(x);

    return ACT_D_exp(logu - exp(logu) - log(x));
}

double pinvexp(double q, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(scale))
	return q + scale;
#endif
    if (!R_FINITE(scale) || scale <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(log(scale) - log(q));

    return ACT_DT_val(exp(-u));
}

double qinvexp(double p, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(scale))
	return p + scale;
#endif
    if (!R_FINITE(scale) || scale <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return -scale / log(ACT_D_Lval(p));
}

double rinvexp(double scale)
{
    if (!R_FINITE(scale) || scale <= 0.0)
        return R_NaN;

    return scale / rexp(1.0);
}

double minvexp(double order, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(scale))
	return order + scale;
#endif
    if (!R_FINITE(scale) ||
        !R_FINITE(order) ||
        scale <= 0.0)
        return R_NaN;

    if (order >= 1.0)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(1.0 - order);
}

double levinvexp(double limit, double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(scale) || ISNAN(order))
	return limit + scale + order;
#endif
    if (!R_FINITE(scale) ||
        R_FINITE(order) ||
        scale <= 0.0)
        return R_NaN;

    if (order >= 1.0)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u = exp(log(scale) - log(limit));

    return R_pow(scale, order) * actuar_gamma_inc(1.0 - order, u)
        + ACT_DLIM__0(limit, order) * (0.5 - exp(-u) + 0.5);
}
