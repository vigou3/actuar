/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-modified negative binomial distribution. See
 *  ../R/ZeroModifiedNegativeBinomial.R for details.
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

/* The negative binomial distribution has p0 = prob^size.
 *
 * Limiting cases:
 *
 * 1. size == 0 is Zero Modified Logarithmic(1 - prob) (according to
 *    the standard parametrization of the logarithmic distribution
 *    used by {d,p,q,r}logarithmic();
 * 2. prob == 1 is mass (1-p0) at x = 1.
 */

double dzmnbinom(double x, double size, double prob, double p0m, int give_log)
{
    /* We compute Pr[X = 0] with dbinom_raw() [as would eventually
     * dnbinom()] to take advantage of all the optimizations for
     * small/large values of 'prob' and 'size' (and also to skip some
     * validity tests).
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0 || !R_FINITE(x)) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* limiting case as size approaches zero is zero modified logarithmic */
    if (size == 0) return dzmlogarithmic(x, 1 - prob, p0m, give_log);

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    return ACT_D_val((1 - p0m) * dnbinom(x, size, prob, /*give_log*/0) / (-expm1(lp0)));
}

double pzmnbinom(double x, double size, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* simple case for all x >= 1 */
    if (p0m == 1) return ACT_DT_1;

    /* limiting case as size approaches zero is zero modified logarithmic */
    if (size == 0) return pzmlogarithmic(x, 1 - prob, p0m, lower_tail, log_p);

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1) return ACT_DT_1;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    return ACT_DT_Cval((1 - p0m) * pnbinom(x, size, prob, /*l._t.*/0, /*log_p*/0)/(-expm1(lp0)));
}

double qzmnbinom(double x, double size, double prob, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(p0m))
	return x + size + prob + p0m;
#endif
    if (prob <= 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as size approaches zero is zero modified logarithmic */
    if (size == 0) return qzmlogarithmic(x, 1 - prob, p0m, lower_tail, log_p);

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1)
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
    x = ACT_DT_1mqIv(x);

    double p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    return qnbinom((1 - p0) * x/(1 - p0m), size, prob, /*l._t.*/0, /*log_p*/0);
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

#define ACT_DIFFMAX_REJECTION 0.6

double rzmnbinom(double size, double prob, double p0m)
{
    if (!R_FINITE(prob) || prob <= 0 || prob > 1 || size < 0 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as size approaches zero is zero modified logarithmic */
    if (size == 0) return rzmlogarithmic(1 - prob, p0m);

    /* limiting case as prob approaches one is mass (1-p0m) at one */
    if (prob == 1) return (unif_rand() <= p0m) ? 0.0 : 1.0;

    double x, p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    /* p0m >= p0: generate from mixture */
    if (p0m >= p0)
	return (unif_rand() * (1 - p0) < (1 - p0m)) ? rnbinom(size, prob) : 0.0;

    /* p0m < p0: choice of algorithm depends on difference p0 - p0m */
    if (p0 - p0m < ACT_DIFFMAX_REJECTION)
    {
	/* rejection method */
	for (;;)
	{
	    x = rnbinom(size, prob);
	    if (x != 0 || /* x == 0 and */ runif(0, p0 * (1 - p0m)) <= (1 - p0) * p0m)
		return x;
	}
    }
    else
    {
    	/* inversion method */
    	return qnbinom(runif((p0 - p0m)/(1 - p0m), 1), size, prob, 1, 0);
    }
}
