/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-modified binomial distribution. See
 *  ../R/ZeroModifiedBinomial.R for details.
 *
 *  Zero-modified distributions are discrete mixtures between a
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
#include "actuar.h"

/* The binomial distribution has p0 = (1 - prob)^size.
 *
 * Limiting cases:
 *
 * 1. size == 0 is mass (1-p0m) at x = 1;
 * 2. prob == 0 is mass (1-p0m) at x = 1.
 */

double dzmbinom(double x, double size, double prob, double p0m, int give_log)
{
    /* We compute Pr[X = 0] with dbinom_raw() [as would eventually
     * dbinom()] to take advantage of all the optimizations for
     * small/large values of 'prob' and 'size' (and also to skip some
     * validity tests).
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob < 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0 || !R_FINITE(x)) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* limiting cases as size -> 1 or prob -> 0 are mass (1-p0m) at one */
    if (size == 1 || prob == 0) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    double lp0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/1);

    return ACT_D_val((1 - p0m) * dbinom(x, size, prob, /*give_log*/0)/(-expm1(lp0)));
}

double pzmbinom(double x, double size, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob < 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* limiting cases as size -> 1 or prob -> 0 are mass (1-p0m) at one */
    if (size == 1 || prob == 0) return ACT_DT_1;

    double lp0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/1);

    return ACT_DT_Cval((1 - p0m) * pbinom(x, size, prob, /*l._t.*/0, /*log_p*/0)/(-expm1(lp0)));
}

double qzmbinom(double x, double size, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob < 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting cases as size -> 1 or prob -> 0 are mass (1-p0m) at one */
    if (size == 1 || prob == 0)
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

    ACT_Q_P01_boundaries(x, 1, size);
    x = ACT_DT_1mqIv(x);

    double p0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/0);

    return qbinom((1 - p0) * x/(1 - p0m), size, prob, /*l._t.*/0, /*log_p*/0);
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

#define ACT_DIFFMAX_REJECTION 0.9

double rzmbinom(double size, double prob, double p0m)
{
    if (!R_FINITE(prob) || prob < 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting cases as size -> 1 or prob -> 0 are mass (1-p0m) at one */
    if (size == 1 || prob == 0) return (unif_rand() <= p0m) ? 0.0 : 1.0;

    double x, p0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/0);

    /* p0m >= p0: generate from mixture */
    if (p0m >= p0)
	return (unif_rand() * (1 - p0) < (1 - p0m)) ? rbinom(size, prob) : 0.0;

    /* p0m < p0: choice of algorithm depends on difference p0 - p0m */
    if (p0 - p0m < ACT_DIFFMAX_REJECTION)
    {
	/* rejection method */
	for (;;)
	{
	    x = rbinom(size, prob);
	    if (x != 0 || /* x == 0 and */ runif(0, p0 * (1 - p0m)) <= (1 - p0) * p0m)
		return x;
	}
    }
    else
    {
    	/* inversion method */
    	return qbinom(runif((p0 - p0m)/(1 - p0m), 1), size, prob, 1, 0);
    }
}
