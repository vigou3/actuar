/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the logarithmic
 *  discrete distribution. See ../R/Logarithmic.R for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dlogarithmic(double x, double p, int give_log)
{
    /*  We work with the probability mass function expressed as
     *
     *  a * p^x / x,  x = 1, 2, ...
     *
     *  with a = -1/log(1 - p).
     */

    /* limiting case as p approaches zero is point mass at one */
    if (p == 0.0)
	return (x == 1.0) ? ACT_D__1 : ACT_D__0;

    if (p < 0.0 || p > 1.0)
        return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 1.0)
        return ACT_D__0;

    x = ACT_forceint(x);

    double a = -1.0/log1p(-p);

    return ACT_D_exp(log(a) + x * log(p) - log(x));
}

double plogarithmic(double x, double p, int lower_tail, int log_p)
{
    /*  We work with the distribution function expressed as
     *
     *  1 - a * pbeta(p, x + 1, 0),
     *
     *  with a = -1/log(1 - p).
     */

    /* limiting case as p approaches zero is point mass at one. */
    if (p == 0.0)
        return (x >= 1.0) ? ACT_DT_1 : ACT_DT_0;

    if (p < 0.0 || p > 1.0)
        return R_NaN;

    if (x <= 1.0)
        return ACT_DT_0;
    if (!R_FINITE(x))
	return ACT_DT_1;

    double a = -1/log1p(-p);

    return ACT_DT_Cval(a * pbeta(p, x + 1, 0, 1, 0));
}

/* double qlogarithmic(double p, double p, int lower_tail, int log_p) */
/* { */
/*     if (p <= 0.0 || p > 1.0) */
/*         return R_NaN; */

/*     ACT_Q_P01_boundaries(p, 1, R_PosInf); */
/*     p =  ACT_D_qIv(p); */

/*     double a; */

/*     a = -1.0/log1p(-p); */

/*     return qbeta(ACT_D_Cval(p)/a, ACT_D_Cval(p) + 1, 0, 1, 0) */
/*     return scale * R_pow(R_pow(ACT_D_Cval(p), -1.0/shape1) - 1.0, 1.0/shape2); */
/* } */

double rlogarithmic(double p)
{
    /*  Implementation of the LS and LK algorithms of:
     *
     *  Kemp, A. W. (1981), Efficient Generation of Logarithmically
     *  Distributed Pseudo-Random Variables, Journal of the Royal
     *  Statistical Society, Series C. Vol. 30, p. 249-253.
     *  URL http://www.jstor.org/stable/2346348
     *
     *  The algorithms are also discussed in chapter 10 of Devroye (1986).
     */

    /* limiting case as p approaches zero is point mass at one. */
    if (p == 0.0)
        return 1.0;

    if (p < 0.0 || p > 1.0)
        return R_NaN;

    /* Automatic selection between the LS and LK algorithms */
    if (p < 0.95)
    {
	double s = -p/log1p(-p);
	double x = 1.0;
	double u = unif_rand();

	while (u > s)
	{
	    u -= s;
	    x += 1.0;
	    s *= p * (x - 1.0)/x;
	}

	return(x);
    }

    /* else (p >= 0.95) */
    double r = log1p(-p);
    double v = unif_rand();

    if (v >= p)     return 1.0;

    double u = unif_rand();
    double q = -expm1(r * u);

    if (v <= (q * q)) return(round(1.0 + log(v)/log(q)));
    if (v <= q)     return(1.0); /* case q^2 < v <= q */
    return(2.0);		   /* case v > q */
}
