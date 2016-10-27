/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, and to simulate random variates for the Poisson-inverse
 *  gaussian distribution. See ../R/PoissonInverseGaussian.R for
 *  details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"
#include "actuar.h"

double dpoisinvgauss(double x, double mu, double phi, int give_log)
{
    /*  We work with the density expressed as
     *
     *  p(x) = sqrt(1/phi) sqrt(1/(pi/2)) exp(1/(phi mu))
     *         * [2 phi (1 + (2 phi mu^2)^(-1))]^(-(x - 0.5)/2)
     *         * besselK(sqrt(2/phi (1 + (2 phi mu^2)^(-1))))
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(phi))
	return x + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;
    ACT_D_nonint_check(x);

    if (!R_FINITE(x) || x < 0.0)
	return ACT_D__0;

    /* limiting case phi = Inf */
    if (!R_FINITE(phi))
	return (x == 0) ? ACT_D__1 : ACT_D__0;

    /* limiting case mu = Inf */
    /* if (!R_FINITE(mu)) */
    /* 	return ACT_D_exp(-(log(phi) + 3 * log(x) + 1/phi/x)/2 - M_LN_SQRT_2PI); */

    /* standard cases */
    double phim = phi * mu, lphi = log(phi);
    double a = 1/(2 * phim * mu), y = x - 0.5;
    double tmp;		      /* log of everything before besselK() */
    double K;		      /* value of the Bessel function */

    tmp = -lphi/2 - M_LN_SQRT_PId2 + 1/phim
	- y * (M_LN2 + lphi + log1p(a))/2 - lgamma(x + 1);

    K = bessel_k(sqrt(2 * (1 + a)/phi), y, /*expo*/1);

    return give_log ? tmp + log(K) : exp(tmp) * K;
}

/*  For ppoisinvgauss(), there does not seem to be algorithms much more
 *  elaborate that successive computations of the probabilities using
 *  the recurrence relationship
 *
 *    p(0) = exp((1 - sqrt(b)))
 *    p(1) = mu p(0) / sqrt(b)
 *
 *  and
 *
 *    p(k) = [a (1 - 1.5/k) p(k - 1) + mu^2/(k (k - 1)) p(k - 2)] / b,
 *
 *  for k = 2, 3, ..., with a = 2 phi mu^2, b = 1 + a.
 */

double ppoisinvgauss(double q, double mu, double phi, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(mu) || ISNAN(phi))
	return q + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    if (q < 0)
        return ACT_DT_0;

    /* limiting case phi = Inf */
    /* if (!R_FINITE(phi)) */
    /* 	return ACT_DT_1; */

    if (!R_FINITE(q))
	return ACT_DT_1;

    /* limiting case mu = Inf */
    /* if (!R_FINITE(mu)) */
    /* 	return pchisq(1/q/phi, 1, !lower_tail, log_p); */

    int k;
    double s, pk, pkm1, pkm2;
    double phim = phi * mu;
    double a = 2 * phim * mu, b = 1 + a;
    double sqr = sqrt(b), mu2 = mu * mu;

    pkm1 = exp((1 - sqr)/phim);	/* p(0) */
    s = pkm1;			/* F(0) */
    if (q == 0) return ACT_DT_val(s);

    pk = (mu/sqr) * pkm1;	/* p(1) */
    s += pk;			/* F(1) */
    if (q == 1) return ACT_DT_val(s);

    for (k = 2; k <= q; k++)
    {
	pkm2 = pkm1;
	pkm1 = pk;
	pk = (a * (1-1.5/k) * pkm1 + mu2/k/(k-1) * pkm2)/b;
	s += pk;
    }

    return ACT_DT_val(s);
}

/* For qpoiinvgauss(), we mostly reuse the code for qnbinom() et al.
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
do_search(double y, double *z, double p, double mu, double phi, double incr)
{
    if(*z >= p) {	/* search to the left */
	for(;;) {
	    if(y == 0 ||
	       (*z = ppoisinvgauss(y - incr, mu, phi, /*l._t.*/1, /*log_p*/0)) < p)
		return y;
	    y = fmax2(0, y - incr);
	}
    }
    else {		/* search to the right */
	for(;;) {
	    y = y + incr;
	    if((*z = ppoisinvgauss(y, mu, phi, /*l._t.*/1, /*log_p*/0)) >= p)
		return y;
	}
    }
}

double qpoisinvgauss(double p, double mu, double phi, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(phi))
	return p + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    /* limiting case phi = Inf */
    /* if (!R_FINITE(phi)) */
    /* 	return 1.0; */

    /* limiting case mu = Inf */
    if (!R_FINITE(mu))
	return 1/phi/qchisq(p, 1, !lower_tail, log_p);

    ACT_Q_P01_boundaries(p, 0, R_PosInf);

    double sigma, sigma2, gamma, z, y;
    double phim = phi * mu;

    /* mu = mu; */
    sigma = mu * sqrt(phim + 1);
    sigma2 = sigma * sigma;
    gamma = (mu + 3*sigma2*(sigma2/mu - 1))/sigma2/sigma;

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

    z = ppoisinvgauss(y, mu, phi, /*lower_tail*/1, /*log_p*/0);

    /* ## From R sources ##
     * fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;

    /* ## From R sources ##
     * If the C-F value is not too large a simple search is OK */
    if (y < 1e5) return do_search(y, &z, p, mu, phi, 1);
    /* ## From R sources ##
     * Otherwise be a bit cleverer in the search */
    {
	double incr = floor(y * 0.001), oldincr;
	do {
	    oldincr = incr;
	    y = do_search(y, &z, p, mu, phi, incr);
	    incr = fmax2(1, floor(incr/100));
	} while(oldincr > 1 && incr > y*1e-15);
	return y;
    }
}

double rpoisinvgauss(double mu, double phi)
{
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    return rpois(rinvgauss(mu, phi));
}
