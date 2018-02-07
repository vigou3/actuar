/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-truncated geometric distribution. See
 *  ../R/ZeroTruncatedGeometric.R for details.
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

/* The geometric distribution has
 *
 *   F(0) = Pr[X = 0] = prob.
 *
 * Limiting case: prob == 1 is point mass at x = 1.
 */

double dztgeom(double x, double prob, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (x == 1) ? ACT_D__1 : ACT_D__0;

    return ACT_D_val(dgeom(x - 1, prob, /*give_log*/0));
}

double pztgeom(double x, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting case as prob approaches one is point mass at one */
    if (prob == 1) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    return ACT_DT_Cval(pgeom(x - 1, prob, /*l._t.*/0, /*log_p*/0));
}

double qztgeom(double x, double prob, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(prob))
	return x + prob;
#endif
    if (prob <= 0 || prob > 1) return R_NaN;

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

    return 1 + qgeom(x, prob, /*l._t.*/1, /*log_p*/0);
}

double rztgeom(double prob)
{
    if (!R_FINITE(prob) || prob <= 0 || prob > 1) return R_NaN;

    /* limiting case as p approaches one is point mass at one */
    if (prob == 1) return 1.0;

    return 1 + rpois(exp_rand() * ((1 - prob) / prob));
}
