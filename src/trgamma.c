/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed gamma distribution. See ../R/TransformedGamma.R
 *  for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dtrgamma(double x, double shape1, double shape2, double scale,
                int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape2 * u^shape1 * e^(-u) / (x * gamma(shape1))
     *
     *  with u = (x/scale)^shape2.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return x + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape1 * shape2 < 1) return R_PosInf;
	if (shape1 * shape2 > 1) return ACT_D__0;
	/* else */
	return give_log ?
	    log(shape2) - log(scale) - lgammafn(shape1) :
	    shape2 / (scale * gammafn(shape1));
    }

    double logu = shape2 * (log(x) - log(scale));

    return ACT_D_exp(log(shape2) + shape1 * logu - exp(logu)
                   - log(x) - lgammafn(shape1));
}

double ptrgamma(double q, double shape1, double shape2, double scale,
                int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return q + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(shape2 * (log(q) - log(scale)));

    return pgamma(u, shape1, 1.0, lower_tail, log_p);
}

double qtrgamma(double p, double shape1, double shape2, double scale,
                int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return p + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * R_pow(qgamma(p, shape1, 1.0, lower_tail, 0),
                         1.0 / shape2);
}

double rtrgamma(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    return scale * R_pow(rgamma(shape1, 1.0), 1.0 / shape2);
}

double mtrgamma(double order, double shape1, double shape2, double scale,
                int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return order + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (order <= -shape1 * shape2)
	return R_PosInf;

    return R_pow(scale, order) * gammafn(shape1 + order/shape2)
        / gammafn(shape1);
}

double levtrgamma(double limit, double shape1, double shape2, double scale,
                  double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= -shape1 * shape2)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u, tmp;

    tmp = shape1 + order / shape2;

    u = exp(shape2 * (log(limit) - log(scale)));

    return R_pow(scale, order) * gammafn(tmp)
        * pgamma(u, tmp, 1.0, 1, 0) / gammafn(shape1)
        + ACT_DLIM__0(limit, order) * pgamma(u, shape1, 1.0, 0, 0);
}
