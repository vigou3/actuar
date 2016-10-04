/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero modified Poisson distribution. See ../R/ZeroModifiedPoisson.R
 *  for details.
 *
 *  Zero modified distributions are discrete mixtures between a
 *  degenerate distribution at zero and the corresponding,
 *  non-modified, distribution. As a mixture, they have density
 *
 *      Pr[Z = x] = [1 - (1 - p0m)/(1 - p0)] 1(x)
 *                  + [(1 - p0m)/(1 - p0)] Pr[X = 0],
 *
 *  where p0 = Pr[X = 0]. The density can also be expressed as
 *  Pr[Z = 0] = p0m and
 *
 *      Pr[Z = x] = (1 - p0m) * Pr[X = x]/(1 - Pr[X = 0]),
 *
 *  for x = 1, 2, ... The distribution function is, for all x,
 *
 *      Pr[Z <= x] = p0m +
 *           (1 - p0m) * (Pr[X <= x] - p0)/(1 - p0)
 *
 *  or, alternatively, the survival function is
 *
 *      Pr[Z > x] = (1 - p0m) * Pr[X > x]/(1 - p0).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

/* The Poisson distribution has p0 = exp(-lambda).
 *
 * Limiting case: lambda == 0 has mass (1 - p0m) at x = 1.
 */

double dzmpois(double x, double lambda, double p0m, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda) || ISNAN(p0m))
	return x + lambda + p0m;
#endif
    if (lambda < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0 || !R_FINITE(x)) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* simple case for all x > 0 */
    if (p0m == 1) return ACT_D__0; /* for all x > 0 */

    /* limiting case as lambda approaches zero is mass (1-p0m) at one */
    if (lambda == 0) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    return ACT_D_exp(dpois(x, lambda, /*give_log*/1)
		     + log1p(-p0m) - ACT_Log1_Exp(-lambda));
}

double pzmpois(double x, double lambda, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda) || ISNAN(p0m))
	return x + lambda + p0m;
#endif
    if (lambda < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* simple case for all x >= 1 */
    if (p0m == 1) return ACT_DT_1;

    /* limiting case as lambda approaches zero is mass (1-p0m) at one */
    if (lambda == 0) return ACT_DT_1;

    return ACT_DT_Cval((1 - p0m) * ppois(x, lambda, /*l._t.*/0, /*log_p*/0) / (-expm1(-lambda)));
}

double qzmpois(double x, double lambda, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda) || ISNAN(p0m))
	return x + lambda + p0m;
#endif
    if (lambda < 0 || !R_FINITE(lambda) || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as lambda approaches zero is mass (1-p0m) at one */
    if (lambda == 0)
    {
	/* simplified ACT_Q_P01_boundaries macro */
	if (log_p)
	{
	    if (x > 0)
		return R_NaN;
	    return (x <= log(p0m)) ? 0.0 : 1.0;
	}
	else /* !log_p */
	{
	    if (x < 0 || x > 1)
		return R_NaN;
	    return (x <= p0m) ? 0.0 : 1.0;
	}
    }

    ACT_Q_P01_boundaries(x, 0, R_PosInf);
    x = ACT_DT_qIv(x);

    double p0c = -expm1(-lambda); /* 1 - Pr[Poisson = 0] */

    return qpois(p0c * (x - p0m)/(1 - p0m) + (1 - p0c), lambda, /*l._t.*/1, /*log_p*/0);
}

/* ALGORITHM FOR GENERATION OF RANDOM VARIATES
 *
 * 1. p0m >= p0: just simulate variates from the discrete mixture.
 *
 * 2. p0m < p0: fastest method depends on the difference p0 - p0m.
 *
 *    2.1 p0 - p0m < ACT_DIFFMAX_REJECTION: rejection method with an
 *        envelope that differs from the target distribution at zero
 *        only. In other words: rejection only at zero.
 *    2.2 p0 - p0m >= ACT_DIFFMAX_REJECTION: inverse method on a
 *        restricted range --- same method as the corresponding zero
 *        truncated distribution.
 *
 * The threshold ACT_DIFFMAX_REJECTION is distribution specific.
 */

#define ACT_DIFFMAX_REJECTION 0.95

double rzmpois(double lambda, double p0m)
{
    if (lambda < 0 || !R_FINITE(lambda) || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as lambda approaches zero is mass (1-p0m) at one */
    if (lambda == 0) return (runif(0, 1) <= p0m) ? 0.0 : 1.0;

    double x, p0 = exp(-lambda);

    /* p0m >= p0: generate from mixture */
    if (p0m >= p0)
	return (runif(0, 1) * (1 - p0) < (1 - p0m)) ? rpois(lambda) : 0.0;

    /* p0m < p0: choice of algorithm depends on difference p0 - p0m */
    if (p0 - p0m < ACT_DIFFMAX_REJECTION)
    {
	/* rejection method */
	for (;;)
	{
	    x = rpois(lambda);
	    if (x != 0 || /* x == 0 and */ runif(0, p0 * (1 - p0m)) <= (1 - p0) * p0m)
		return x;
	}
    }
    else
    {
	/* inversion method */
	return qpois(runif((p0 - p0m)/(1 - p0m), 1), lambda, 1, 0);
    }
}
