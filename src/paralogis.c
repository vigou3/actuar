/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the paralogistic distribution. See ../R/Paralogistic.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dparalogis(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape^2 * u^shape * (1 - u) / x
     *
     *  with u = 1/(1 + v), v = (x/scale)^shape.
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

    double tmp, logu, log1mu;

    tmp = shape * (log(x) - log(scale));
    logu = - log1pexp(tmp);
    log1mu = - log1pexp(-tmp);

    return ACT_D_exp(2.0 * log(shape) + shape * logu + log1mu - log(x));
}

double pparalogis(double q, double shape, double scale, int lower_tail,
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
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    double u = exp(-log1pexp(shape * (log(q) - log(scale))));

    return ACT_DT_Cval(R_pow(u, shape));
}

double qparalogis(double p, double shape, double scale, int lower_tail,
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
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    double tmp = 1.0 / shape;

    return scale * R_pow(R_pow(ACT_D_Cval(p), -tmp) - 1.0, tmp);
}

double rparalogis(double shape, double scale)
{
    double tmp;

    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    tmp = 1.0 / shape;

    return scale * R_pow(R_pow(unif_rand(), -tmp) - 1.0, tmp);
}

double mparalogis(double order, double shape, double scale, int give_log)
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

    if (order <= -shape ||
        order >= shape * shape)
        return R_PosInf;

    double tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(shape - tmp)
        / gammafn(shape);
}

double levparalogis(double limit, double shape, double scale, double order,
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

    if (order <= -shape)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double tmp = order / shape;
    double u = exp(-log1pexp(shape * (log(limit) - log(scale))));

    return R_pow(scale, order)
	* betaint_raw(0.5 - u + 0.5, 1.0 + tmp, shape - tmp)
	/ gammafn(shape)
        + ACT_DLIM__0(limit, order) * R_pow(u, shape);
}
