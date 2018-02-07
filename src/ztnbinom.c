/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-truncated negative binomial distribution. See
 *  ../R/ZeroTruncatedNegativeBinomial.R for details.
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
#include "actuar.h"

/* The negative binomial distribution has
 *
 *   F(0) = Pr[X = 0] = prob^size.
 *
 * Limiting cases:
 *
 * 1. size == 0 is Logarithmic(1 - prob) (according to the standard
 *    parametrization of the logarithmic distribution used by
 *    {d,p,q,r}logarithmic();
 * 2. prob == 1 is point mass at x = 1.
 */

double dztnbinom(double x, double size, double prob, int give_log)
{
    /* We compute Pr[X = 0] with dbinom_raw() [as would eventually
     * dnbinom()] to take advantage of all the optimizations for
     * small/large values of 'prob' and 'size' (and also to skip some
     * validity tests).
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob <= 0 || prob > 1 || size < 0) return R_NaN;

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting case as size approaches zero is logarithmic */
    if (size == 0) return dlogarithmic(x, 1 - prob, give_log);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (x == 1) ? ACT_D__1 : ACT_D__0;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    return ACT_D_val(dnbinom(x, size, prob, /*give_log*/0)/(-expm1(lp0)));
}

double pztnbinom(double x, double size, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob <= 0 || prob > 1 || size < 0) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting case as size approaches zero is logarithmic */
    if (size == 0) return plogarithmic(x, 1 - prob, lower_tail, log_p);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    return ACT_DT_Cval(pnbinom(x, size, prob, /*l._t.*/0, /*log_p*/0)/(-expm1(lp0)));
}

double qztnbinom(double x, double size, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(size) || ISNAN(prob))
	return x + size + prob;
#endif
    if (prob <= 0 || prob > 1 || size < 0) return R_NaN;

    /* limiting case as size approaches zero is logarithmic */
    if (size == 0) return qlogarithmic(x, 1 - prob, lower_tail, log_p);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1)
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

    ACT_Q_P01_boundaries(x, 1, R_PosInf);
    x = ACT_DT_qIv(x);

    double p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    return qnbinom(p0 + (1 - p0) * x, size, prob, /*l._t.*/1, /*log_p*/0);
}

double rztnbinom(double size, double prob)
{
    if (!R_FINITE(prob) || prob <= 0 || prob > 1 || size < 0) return R_NaN;

    /* limiting case as size approaches zero is logarithmic */
    if (size == 0) return rlogarithmic(1 - prob);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return 1.0;

    double p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    return qnbinom(runif(p0, 1), size, prob, /*l._t.*/1, /*log_p*/0);
}
