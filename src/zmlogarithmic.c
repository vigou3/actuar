/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the zero modified
 *  logarithmic discrete distribution. See
 *  ../R/ZeroModifiedLogarithmic.R for details.
 *
 *  The zero modified logarithmic is a discrete mixtures between a
 *  degenerate distribution at zero and a logarithmic distribution.
 *  The density is
 *
 *      Pr[Z = x] = p0m 1(x) + (1 - p0m) Pr[X = x]
 *
 *  or, alternatively, Pr[Z = 0] = p0m and
 *
 *      Pr[Z = x] = (1 - p0m) Pr[X = x],
 *
 *  for x = 1, 2, ... The distribution function is, for all x,
 *
 *      Pr[Z <= x] = p0m + (1 - p0m) Pr[X <= x]
 *
 *  and the survival function is
 *
 *      Pr[Z > x] = (1 - p0m) * Pr[X > x].
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dzmlogarithmic(double x, double p, double p0m, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p) || ISNAN(p0m))
	return x + p + p0m;
#endif
    if (p < 0 || p >= 1 || p0m < 0 || p0m > 1) return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 0) return ACT_D__0;
    if (x == 0) return ACT_D_val(p0m);
    /* NOTE: from now on x > 0 */

    /* limiting case as p approaches zero is mass (1-p0m) at one */
    if (p == 0) return (x == 1) ? ACT_D_Clog(p0m) : ACT_D__0;

    x = ACT_forceint(x);

    double a = -1.0/log1p(-p);

    return ACT_D_exp(log(a) + x * log(p) + log1p(-p0m) - log(x));
}

double pzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p) || ISNAN(p0m))
	return x + p + p0m;
#endif
    if (p < 0 || p >= 1 || p0m < 0 || p0m > 1) return R_NaN;

    if (x < 0) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;
    if (x < 1) return ACT_DT_val(p0m);
    /* NOTE: from now on x >= 1 */

    /* simple case for all x >= 1 */
    if (p0m == 1) return ACT_DT_1;

    /* limiting case as p approaches zero is mass (1-p0m) at one. */
    if (p == 0) return ACT_DT_1;

    return ACT_DT_Cval((1 - p0m) * plogarithmic(x, p, /*l._t.*/0, /*log_p*/0));
}

double qzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p) || ISNAN(p0m))
	return x + p + p0m;
#endif
    if (p < 0 || p >= 1 || p0m < 0 || p0m > 1) return R_NaN;

    /* limiting case as p approaches zero is mass (1-p0m) at one */
    if (p == 0)
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

    ACT_Q_P01_boundaries(x, 1.0, R_PosInf);
    x = ACT_DT_qIv(x);

    /* avoid rounding errors below if x was given in log form */
    if (log_p)
	p0m = exp(log(p0m));

    return (x <= p0m) ? 0.0 : qlogarithmic((x - p0m)/(1 - p0m), p, /*l._t.*/1, /*log_p*/0);
}

/* ALGORITHM FOR GENERATION OF RANDOM VARIATES
 *
 * Just simulate variates from the discrete mixture.
 *
 */
double rzmlogarithmic(double p, double p0m)
{
    if (p < 0 || p >= 1 || p0m < 0 || p0m > 1) return R_NaN;

    return (runif(0, 1) < p0m) ? 0.0 : rlogarithmic(p);
}
