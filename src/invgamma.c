/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Inverse Gamma distribution. See ../R/InverseGamma.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dinvgamma(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  u^shape * e^(-u) / (x * gamma(shape))
     *
     *  with u = scale/x.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
	return x + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    /* handle also x == 0 here */
    if (!R_FINITE(x) || x <= 0.0)
        return ACT_D__0;

    double logu = log(scale) - log(x);

    return ACT_D_exp(shape * logu - exp(logu) - log(x) - lgammafn(shape));
}

double pinvgamma(double q, double shape, double scale, int lower_tail,
                 int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape) || ISNAN(scale))
	return q + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(log(scale) - log(q));

    return pgamma(u, shape, 1.0, !lower_tail, log_p);
}

double qinvgamma(double p, double shape, double scale, int lower_tail,
                 int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape) || ISNAN(scale))
	return p + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale / qgamma(p, shape, 1.0, !lower_tail, 0);
}

double rinvgamma(double shape, double scale)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;;

    return scale / rgamma(shape, 1.0);
}

double minvgamma(double order, double shape, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape) || ISNAN(scale))
	return order + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order >= shape)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(shape - order) / gammafn(shape);
}

double levinvgamma(double limit, double shape, double scale, double order,
                   int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape) || ISNAN(scale) || ISNAN(order))
	return limit + shape + scale + order;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order >= shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u = exp(log(scale) - log(limit));

    return R_pow(scale, order) * actuar_gamma_inc(shape - order, u) / gammafn(shape)
        + ACT_DLIM__0(limit, order) * pgamma(u, shape, 1.0, 1, 0);
}

double mgfinvgamma(double t, double shape, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(t) || ISNAN(shape) || ISNAN(scale))
	return t + shape + scale;
#endif
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0 ||
        t > 0.0 )
        return R_NaN;

    if (t == 0.0)
        return ACT_D__1;

    /* rescale and change sign */
    t = -scale * t;

    return ACT_D_exp(M_LN2 + 0.5 * shape * log(t) +
		     log(bessel_k(sqrt(4 * t), shape, 1)) -
		     lgammafn(shape));
}
