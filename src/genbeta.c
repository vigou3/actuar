/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the generalized beta distribution. See ../R/GeneralizedBeta.R
 *  for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dgenbeta(double x, double shape1, double shape2, double shape3,
                double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape3 * u^shape1 * (1 - u)^(shape2 - 1) / (x * beta(shape1, shape2))
     *
     *  with u = (x/scale)^shape3.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return x + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (x < 0.0 || x > scale)
        return ACT_D__0;
    if (x == 0.0)
    {
        double tmp = shape1 * shape3;
        if (tmp > 1) return(ACT_D__0);
        if (tmp < 1) return(R_PosInf);
        /* tmp == 1 : */ return(ACT_D_val(shape3/beta(shape1, shape2)));
    }
    if (x == scale)
    {
        if (shape2 > 1) return(ACT_D__0);
        if (shape2 < 1) return(R_PosInf);
        /* shape2 == 1 : */ return(ACT_D_val(shape1 * shape3));
    }

    double logu, log1mu;

    logu = shape3 * (log(x) - log(scale));
    log1mu = log1p(-exp(logu));

    return ACT_D_exp(log(shape3) + shape1 * logu + (shape2 - 1.0) * log1mu
                   - log(x) - lbeta(shape1, shape2));
}

double pgenbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return q + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;
    if (q >= scale)
        return ACT_DT_1;

    double u = exp(shape3 * (log(q) - log(scale)));

    return pbeta(u, shape1, shape2, lower_tail, log_p);
}

double qgenbeta(double p, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return p + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, scale);
    p = ACT_D_qIv(p);

    return scale * R_pow(qbeta(p, shape1, shape2, lower_tail, 0),
                         1.0 / shape3);
}

double rgenbeta(double shape1, double shape2, double shape3, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * R_pow(rbeta(shape1, shape2), 1.0 / shape3);
}

double mgenbeta(double order, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale))
	return order + shape1 + shape2 + shape3 + scale;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order  <= - shape1 * shape3)
	return R_PosInf;

    double tmp = order / shape3;

    return R_pow(scale, order) * beta(shape1 + tmp, shape2)
        / beta(shape1, shape2);
}

double levgenbeta(double limit, double shape1, double shape2, double shape3,
                 double scale, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(shape1) || ISNAN(shape2) || ISNAN(shape3) || ISNAN(scale) || ISNAN(order))
	return limit + shape1 + shape2 + shape3 + scale + order;
#endif
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order <= - shape1 * shape3)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    double u, tmp;

    tmp = order / shape3;

    u = exp(shape3 * (log(limit) - log(scale)));

    return R_pow(scale, order) * beta(shape1 + tmp, shape2)
        / beta(shape1, shape2) * pbeta(u, shape1 + tmp, shape2, 1, 0)
        + ACT_DLIM__0(limit, order) * pbeta(u, shape1, shape2, 0, 0);
}
