/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero truncated negative binomial distribution. See
 *  ../R/ZeroTruncatedNegativeBinomial.R for details.
 *
 *  Zero truncated distributions have density
 *
 *      Pr[Z = x] = Pr[X = x]/(1 - Pr[X = 0]),
 *
 *  and distribution function
 *
 *      Pr[Z <= x] = (Pr[X <= x] - Pr[X = 0])/(1 - Pr[X = 0])
 *
 *  or, alternatively, survival function
 *
 *      Pr[Z > x] = (1 - Pr[X > x])/(1 - Pr[X = 0]).
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

/* The Zero truncated negative binomial distribution has
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
    /* We compute Pr[X = 0] dbinom_raw() [as would eventually dnbinom()]
     * to take advantage of all the optimizations for small/large values of
     * 'prob' and 'size' (and also to skip some validity tests).
     */

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting case as size approches zero is logarithmic */
    if (size == 0) return dlogarithmic(x, 1 - prob, give_log);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (x == 1) ? ACT_D__1 : ACT_D__0;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    /* limiting case as prob approches 1 handled automatically */
    return ACT_D_val(dnbinom(x, size, prob, /*give_log*/0)/(-expm1m(lp0)));
}

double pztnbinom(double q, double size, double prob, int lower_tail, int log_p)
{
    if (q < 1) return ACT_DT_0;
    if (!R_FINITE(q)) return ACT_DT_1;

    /* limiting case as size approches zero is logarithmic */
    if (size == 0) return plogarithmic(q, 1 - prob, lower_tail, log_p);

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (q >= 1) ? ACT_D__1 : ACT_D__0;

    double lp0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/1);

    return ACT_DT_Cval(pnbinom(q, size, prob, /*l._t.*/0, /*log_p*/0)/(-expm1(lp0)));
}

double qztnbinom(double p, double size, double prob, int lower_tail, int log_p)
{
    /* limiting case as size approches zero is logarithmic */
    if (size == 0) return qlogarithmic(p, 1 - prob, lower_tail, log_p);

    /* limiting case as p approaches one is point mass at one */
    if (prob == 1)
    {
	/* simplified ACT_Q_P01_boundaries macro */
	if (log_p)
	{
	    if (p > 0)
		return R_NaN;
	    return 1.0;
	}
	else /* !log_p */
	{
	    if (p < 0 || p > 1)
		return R_NaN;
	    return 1.0;
	}
    }

    ACT_Q_P01_boundaries(p, 1, R_PosInf);
    p = ACT_D_qIv(p);

    double p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    return qnbinom(p + p0 * (0.5 - p + 0.5), size, prob, /*l._t.*/1, /*log_p*/0);
}

double rnbinom(double size, double prob)
{
    /* limiting case as size approches zero is logarithmic */
    if (size == 0) return rlogarithmic(1 - prob);

    /* limiting case as p approaches one is point mass at one */
    if (prob == 1) return 1.0;

    double p0 = dbinom_raw(size, size, prob, 1 - prob, /*give_log*/0);

    return qnbinom(runif(p0, 1), size, prob, /*l._t.*/1, /*log_p*/0);
}
