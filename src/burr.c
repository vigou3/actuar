/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Burr distribution. See ../R/Burr.R for details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dburr(double x, double shape1, double shape2, double scale,
             int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape1 * shape2 * u^shape1 * (1 - u) / x
     *
     *  with u = 1/(1 + v), v = (x/scale)^shape2.
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
	if (shape2 < 1) return R_PosInf;
	if (shape2 > 1) return ACT_D__0;
	/* else */
	return ACT_D_val(shape1 / scale);
    }

    double tmp, logu, log1mu;

    tmp = shape2 * (log(x) - log(scale));
    logu = - log1pexp(tmp);
    log1mu = - log1pexp(-tmp);

    return ACT_D_exp(log(shape1) + log(shape2) + shape1 * logu
                   + log1mu - log(x));
}

double pburr(double q, double shape1, double shape2, double scale,
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

    double u = exp(-log1pexp(shape2 * (log(q) - log(scale))));

    return ACT_DT_Cval(R_pow(u, shape1));
}

double qburr(double p, double shape1, double shape2, double scale,
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
    p =  ACT_D_qIv(p);

    return scale * R_pow(R_pow(ACT_D_Cval(p), -1.0/shape1) - 1.0, 1.0/shape2);
}

double rburr(double shape1, double shape2, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    return scale * R_pow(R_pow(unif_rand(), -1.0/shape1) - 1.0, 1.0/shape2);
}

double mburr(double order, double shape1, double shape2, double scale,
             int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale))
	return order + shape1 + shape2 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order  <= -shape2 ||
        order  >= shape1 * shape2)
	return R_PosInf;

    double tmp = order / shape2;

    return R_pow(scale, order) * gammafn(1.0 + tmp) * gammafn(shape1 - tmp)
        / gammafn(shape1);
}

double levburr(double limit, double shape1, double shape2, double scale,
               double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= -shape2)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u = exp(-log1pexp(shape2 * (log(limit) - log(scale))));
    double tmp = order / shape2;

    return R_pow(scale, order)
	* betaint_raw(0.5 - u + 0.5, 1.0 + tmp, shape1 - tmp)
	/ gammafn(shape1)
	+ ACT_DLIM__0(limit, order) * R_pow(u, shape1);
}
