/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the Gumbel distribution. See ../R/Gumbel.R for
 *  details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dgumbel(double x, double alpha, double beta, int give_log)
{
    /*  We work with the density expressed as
     *
     *  e^(-u) * e^(-e^(-u)) / beta
     *
     *  with u = (x - alpha)/beta.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
	return x + alpha + beta;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(beta) ||
        beta <= 0.0)
        return R_NaN;;

    double u = (x - alpha)/beta;

    return ACT_D_exp(-(u + exp(-u) + log(beta)));
}

double pgumbel(double q, double alpha, double beta, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(alpha) || ISNAN(beta))
	return q + alpha + beta;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(beta) ||
        beta <= 0.0)
        return R_NaN;;

    double u = (q - alpha)/beta;

    return ACT_DT_val(exp(-exp(-u)));
}

double qgumbel(double p, double alpha, double beta, int lower_tail,
                   int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(alpha) || ISNAN(beta))
	return p + alpha + beta;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(beta) ||
        beta <= 0.0)
        return R_NaN;;

    ACT_Q_P01_boundaries(p, R_NegInf, R_PosInf);
    p = ACT_DT_qIv(p);

    return alpha - beta * log(-log(p));
}

double rgumbel(double alpha, double beta)
{
    if (!R_FINITE(alpha) ||
        !R_FINITE(beta) ||
        beta <= 0.0)
        return R_NaN;;

    return alpha - beta * log(exp_rand());
}

#define EULER_CNST 0.577215664901532860606512090082

double mgumbel(double order, double alpha, double beta, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(alpha) || ISNAN(beta))
	return order + alpha + beta;
#endif
    if (!R_FINITE(alpha) ||
        !R_FINITE(beta) ||
        !R_FINITE(order) ||
        beta <= 0.0 ||
        order <= 0.0 ||
	order > 2.0)
        return R_NaN;

    if (order == 1.0)
	return alpha + EULER_CNST * beta;
    if (order == 2.0)
	return R_pow_di(M_PI * beta, 2)/6 + R_pow_di(alpha + EULER_CNST * beta, 2);

    return R_NaN;		/* order != 1 or 2 */
}
