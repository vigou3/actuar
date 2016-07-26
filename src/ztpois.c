/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute probability function, cumulative distribution
 *  and quantile functions, and to simulate random variates for the
 *  zero truncated Poisson distribution. See
 *  ../R/ZeroTruncatedPoisson.R for details.
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

/* The Poisson distribution has F(0) = Pr[X = 0] = exp(-lambda) */

double dztpois(double x, double lambda, int give_log)
{
    if (x < 1.) return ACT_D__0;

    return ACT_D_exp(dpois(x, lambda, 1) - ACT_Log1_Exp(-lambda));
}

double pztpois(double q, double lambda, int lower_tail, int log_p)
{
    if (q < 1.) return ACT_DT_0;

    return ACT_DT_Cval(ppois(q, lambda, 0, 0)/(-expm1(-lambda)));
}

double qztpois(double p, double lambda, int lower_tail, int log_p)
{
    ACT_Q_P01_boundaries(p, 1, R_PosInf);
    p = ACT_D_qIv(p);

    return qpois(p + exp(-lambda) * (0.5 - p + 0.5), lambda, 1, 0);
}

double rztpois(double lambda)
{
    if (lambda <= 0.) return 1.;

    return qpois(runif(exp(-lambda), 1), lambda, 1, 0);
}
