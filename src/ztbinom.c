/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-truncated binomial distribution. See
 *  ../R/ZeroTruncatedBinomial.R for details.
 *
 *  Zero-truncated distributions have density
 *
 *      Pr[Z = x] = Pr[X = x]/(1 - Pr[X = 0]),
 *
 *  and distribution function
 *
 *      Pr[Z <= x] = (Pr[X <= x] - Pr[X = 0])/(1 - Pr[X = 0])
 *
 *  or, alternatively, survival function
 *
 *      Pr[Z > x] = Pr[X > x]/(1 - Pr[X = 0]).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

/* The binomial distribution has
 *
 *   F(0) = Pr[X = 0] = (1 - prob)^size.
 *
 * Support is x = 1, ..., size.
 *
 * Limiting cases:
 *
 * 1. size == 1 is point mass at x = 1;
 * 2. prob == 0 is point mass at x = 1.
 */

double dztbinom(double x, double size, double prob, int give_log)
{
    /* We compute Pr[X = 0] with dbinom_raw() [as would eventually
     * dbinom()] to take advantage of all the optimizations for
     * small/large values of 'prob' and 'size' (and also to skip some
     * validity tests).
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob < 0 || prob > 1 || size < 1) return R_NaN;

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting cases as size -> 1 or prob -> 0 are point mass at one */
    if (size == 1 || prob == 0) return (x == 1) ? ACT_D__1 : ACT_D__0;

    double lp0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/1);

    return ACT_D_val(dbinom(x, size, prob, /*give_log*/0)/(-expm1(lp0)));
}

double pztbinom(double x, double size, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob < 0 || prob > 1 || size < 1) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting cases as size -> 1 or prob -> 0 are point mass at one */
    if (size == 1 || prob == 0) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    double lp0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/1);

    return ACT_DT_Cval(pbinom(x, size, prob, /*l._t.*/0, /*log_p*/0)/(-expm1(lp0)));
}

double qztbinom(double x, double size, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob < 0 || prob > 1 || size < 1) return R_NaN;

    /* limiting cases as size -> 1 or prob -> 0 are point mass at one */
    if (size == 1 || prob == 0)
    {
	/* simplified ACT_Q_P01_boundaries macro */
	if (log_p)
	{
	    if (x > 0)
		return R_NaN;
	    return 1.0;
	}
	else /* !log_p */
	{
	    if (x < 0 || x > 1)
		return R_NaN;
	    return 1.0;
	}
    }

    ACT_Q_P01_boundaries(x, 1, size);
    x = ACT_DT_qIv(x);

    double p0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/0);

    return qbinom(p0 + (1 - p0) * x, size, prob, /*l._t.*/1, /*log_p*/0);
}

double rztbinom(double size, double prob)
{
    if (!R_FINITE(prob) || prob < 0 || prob > 1 || size < 0) return R_NaN;

    /* limiting cases as size -> 1 or prob -> 0 are point mass at one */
    if (size == 1 || prob == 0) return 1.0;

    double p0 = dbinom_raw(0, size, prob, 1 - prob, /*give_log*/0);

    return qbinom(runif(p0, 1), size, prob, /*l._t.*/1, /*log_p*/0);
}
