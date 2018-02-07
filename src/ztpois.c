/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero-truncated Poisson distribution. See
 *  ../R/ZeroTruncatedPoisson.R for details.
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

/* The Poisson distribution has
 *
 *   F(0) = Pr[X = 0] = exp(-lambda).
 *
 * Limiting case: lambda == 0 is point mass at x = 1.
 */

double dztpois(double x, double lambda, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda))
	return x + lambda;
#endif
    if (lambda < 0) return R_NaN;

    if (x < 1 || !R_FINITE(x)) return ACT_D__0;

    /* limiting case as lambda approaches zero is point mass at one */
    if (lambda == 0) return (x == 1) ? ACT_D__1 : ACT_D__0;

    return ACT_D_exp(dpois(x, lambda, /*give_log*/1) - ACT_Log1_Exp(-lambda));
}

double pztpois(double x, double lambda, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda))
	return x + lambda;
#endif
    if (lambda < 0) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting case as lambda approaches zero is point mass at one */
    if (lambda == 0) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    return ACT_DT_Cval(ppois(x, lambda, /*l._t.*/0, /*log_p*/0)/(-expm1(-lambda)));
}

double qztpois(double x, double lambda, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(lambda))
	return x + lambda;
#endif
    if (lambda < 0 || !R_FINITE(lambda)) return R_NaN;

    /* limiting case as lambda approaches zero is point mass at one */
    if (lambda == 0)
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

    double p0c = -expm1(-lambda);

    return qpois(p0c * x + (1 - p0c), lambda, /*l._t.*/1, /*log_p*/0);
}

double rztpois(double lambda)
{
    if (lambda < 0 || !R_FINITE(lambda)) return R_NaN;

    /* limiting case as lambda approaches zero is point mass at one */
    if (lambda == 0) return 1.0;

    return qpois(runif(exp(-lambda), 1), lambda, 1, 0);
}
