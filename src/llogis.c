/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the loglogistic distribution. See ../R/Loglogistic.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dllogis(double x, double shape, double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape * u * (1 - u) / x
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape.
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
    logu = - log1pexp(-tmp);
    log1mu = - log1pexp(tmp);

    return ACT_D_exp(log(shape) + logu + log1mu - log(x));
}

double pllogis(double q, double shape, double scale, int lower_tail, int log_p)
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

    double u = exp(-log1pexp(shape * (log(scale) - log(q))));

    return ACT_DT_val(u);
}

double qllogis(double p, double shape, double scale, int lower_tail, int log_p)
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

    return scale * R_pow(1.0 / ACT_D_Cval(p) - 1.0, 1.0/shape);
}

double rllogis(double shape, double scale)
{
    if (!R_FINITE(shape) ||
        !R_FINITE(scale) ||
        shape <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * R_pow(1.0 / unif_rand() - 1.0, 1.0 / shape);
}

double mllogis(double order, double shape, double scale, int give_log)
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
        order >= shape)
	return R_PosInf;

    double tmp = order / shape;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(1.0 - tmp);
}

double levllogis(double limit, double shape, double scale, double order,
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
	return R_NaN;;

    if (order <= -shape)
        return R_PosInf;

    if (limit <= 0.0)
        return 0;

    double u = exp(-log1pexp(shape * (log(scale) - log(limit))));
    double tmp = order / shape;

    return R_pow(scale, order)
	* betaint_raw(u, 1.0 + tmp, 1.0 - tmp)
        + ACT_DLIM__0(limit, order) * (0.5 - u + 0.5);
}
