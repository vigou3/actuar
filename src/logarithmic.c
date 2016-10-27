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

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p))
	return x + p;
#endif
    if (p < 0 || p >= 1) return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 1) return ACT_D__0;

    /* limiting case as p approaches zero is point mass at one */
    if (p == 0) return (x == 1) ? ACT_D__1 : ACT_D__0;

    x = ACT_forceint(x);

    double a = -1.0/log1p(-p);

    return ACT_D_exp(log(a) + x * log(p) - log(x));
}

/*  For plogarithmic(), there does not seem to be algorithms much more
 *  elaborate that successive computations of the probabilities using
 *  the recurrence relationship
 *
 *  P[X = x + 1] = p * x * Pr[X = x] / (x + 1), x = 2, 3, ...
 *
 *  with Pr[X = 1] = -p/log(1 - p). This is what is done here.
 */

double plogarithmic(double x, double p, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p))
	return x + p;
#endif
    if (p < 0 || p >= 1) return R_NaN;

    if (x < 1) return ACT_DT_0;
    if (!R_FINITE(x)) return ACT_DT_1;

    /* limiting case as p approaches zero is point mass at one. */
    if (p == 0) return (x >= 1) ? ACT_DT_1 : ACT_DT_0;

    int k;
    double s, pk;

    pk = -p/log1p(-p);		      /* Pr[X = 1] */
    s = pk;

    if (x == 1) return ACT_DT_val(s); /* simple case */

    for (k = 1; k < x; k++)
    {
	pk *= p * k/(k + 1.0);
	s += pk;
    }

    return ACT_DT_val(s);
}

/* For qlogarithmic(), we mostly reuse the code for qnbinom() et al.
 * in the R sources. From src/nmath/qnbinom.c:
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 */

static double
do_search(double y, double *z, double x, double pr, double incr)
{
    if(*z >= x) {	/* search to the left */
	for(;;) {
	    if(y == 0 ||
	       (*z = plogarithmic(y - incr, pr, /*l._t.*/1, /*log_p*/0)) < x)
		return y;
	    y = fmax2(0, y - incr);
	}
    }
    else {		/* search to the right */
	for(;;) {
	    y = y + incr;
	    if((*z = plogarithmic(y, pr, /*l._t.*/1, /*log_p*/0)) >= x)
		return y;
	}
    }
}

double qlogarithmic(double x, double p, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(p))
	return x + p;
#endif
    if (p < 0 || p >= 1) return R_NaN;

    /* limiting case as p approaches zero is point mass at one */
    if (p == 0)
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

    ACT_Q_P01_boundaries(x, 1.0, R_PosInf);

    double a, P, Q, mu, sigma, gamma, z, y;

    a = -1.0/log1p(-p);
    P = a * p;
    Q = 1.0/(0.5 - p + 0.5);
    mu = P * Q;
    sigma = sqrt(P * (1 - P)) * Q;
    gamma = (P * (1 + p - 3*P + 2*P*P) * R_pow_di(Q, 3))/R_pow_di(sigma, 3);

    /* ## From R sources ##
     * Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if (!lower_tail || log_p)
    {
	p = ACT_DT_qIv(p); /* need check again (cancellation!): */
	if (p == ACT_DT_0) return 0;
	if (p == ACT_DT_1) return R_PosInf;
    }
    /* ## From R sources ##
     * temporary hack --- FIXME --- */
    if (p + 1.01 * DBL_EPSILON >= 1.0) return R_PosInf;

    /* ## From R sources ##
     * y := approx.value (Cornish-Fisher expansion) :  */
    z = qnorm(p, 0.0, 1.0, /*lower_tail*/1, /*log_p*/0);
    y = ACT_forceint(mu + sigma * (z + gamma * (z*z - 1)/6));

    z = plogarithmic(y, p, /*lower_tail*/1, /*log_p*/0);

    /* ## From R sources ##
     * fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    /* ## From R sources ##
     * If the C-F value is not too large a simple search is OK */
    if (y < 1e5) return do_search(y, &z, x, p, 1);
    /* ## From R sources ##
     * Otherwise be a bit cleverer in the search */
    {
	double incr = floor(y * 0.001), oldincr;
	do {
	    oldincr = incr;
	    y = do_search(y, &z, x, p, incr);
	    incr = fmax2(1, floor(incr/100));
	} while(oldincr > 1 && incr > y*1e-15);
	return y;
    }
}

/*  rlogarithmic() is an implementation with automatic selection of
 *  the LS and LK algorithms of:
 *
 *  Kemp, A. W. (1981), Efficient Generation of Logarithmically
 *  Distributed Pseudo-Random Variables, Journal of the Royal
 *  Statistical Society, Series C. Vol. 30, p. 249-253.
 *  URL http://www.jstor.org/stable/2346348
 *
 *  The algorithms are also discussed in chapter 10 of Devroye (1986).
 */

double rlogarithmic(double p)
{
    if (p < 0 || p > 1) return R_NaN;

    /* limiting case as p approaches zero is point mass at one. */
    if (p == 0) return 1.0;

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
    {
	double r = log1p(-p);
	double v = unif_rand();

	if (v >= p)       return 1.0;

	double u = unif_rand();
	double q = -expm1(r * u);

	if (v <= (q * q)) return(round(1.0 + log(v)/log(q)));
	if (v <= q)       return(1.0); /* case q^2 < v <= q */
	return(2.0);		       /* case v > q */
    }
}
