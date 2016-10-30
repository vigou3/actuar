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
     *  p(x) = sqrt(1/phi) sqrt(1/(pi/2)) exp(1/(phi mu))/x!
     *       * [sqrt(2 phi (1 + (2 phi mu^2)^(-1)))]^(-(x - 0.5))
     *       * bessel_k(sqrt(2/phi (1 + (2 phi mu^2)^(-1))), x - 0.5)
     *
     *  or, is essence,
     *
     *  p(x) = A exp(1/(phi mu))/x! B^(-y) bessel_k(B/phi, y)
     *
     *  The limiting case mu = Inf is handled "automatically" with
     *  terms going to zero when mu is Inf. Specific code not worth it
     *  since the function should rarely be evaluated with mu = Inf in
     *  practice.
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

    /* standard cases (and limiting case mu = Inf) */
    double phim = phi * mu, lphi = log(phi);
    double a = 1/(2 * phim * mu), y = x - 0.5;
    double logA, logB;	      /* log of big constants */
    double lpx;		      /* log of everything before bessel_k() */
    double K;		      /* value of the Bessel function */

    logA = -lphi/2 - M_LN_SQRT_PId2 + 1/phim;
    logB = (M_LN2 + lphi + log1p(a))/2;

    lpx = logA - y * logB - lgamma1p(x);
    K = bessel_k(exp(logB - lphi), y, /*expo*/1);

    return give_log ? lpx + log(K) : exp(lpx) * K;
}

/*  For ppoisinvgauss(), there does not seem to be algorithms much
 *  more elaborate that successive computations of the probabilities.
 *  Performance wise, the explicit formula used in dpoisinvgauss() is
 *  on par (and slightly faster in certain cases) with the recursive
 *  formulas and we get "for free" the optimizations in bessel_k().
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
    if (!R_FINITE(phi))
    	return ACT_DT_1;

    if (!R_FINITE(q))
	return ACT_DT_1;

    int x;
    double phim = phi * mu, lphi = log(phi);
    double a = 1/(2 * phim * mu);
    double logA, logB, C;
    double y, s = 0;

    logA = -lphi/2 - M_LN_SQRT_PId2 + 1/phim;
    logB = (M_LN2 + lphi + log1p(a))/2;
    C = exp(logB - lphi);

    for (x = 0; x <= q; x++)
    {
	y = x - 0.5;
	s += exp(logA - y * logB - lgamma1p(x)) * bessel_k(C, y, /*expo*/1);
    }

    return ACT_D_val(s);
}

/*  For qpoiinvgauss(), we mostly reuse the code for qnbinom() et al.
 *  in the R sources. From src/nmath/qnbinom.c:
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 *
 *  For the limiting case mu = Inf (that has no finite moments), we
 *  use instead the quantile of an inverse chi-square distribution as
 *  starting point.
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
    if (!R_FINITE(phi))
    	return 0.0;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);

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

    double z, y;

    z = qnorm(p, 0.0, 1.0, /*lower_tail*/1, /*log_p*/0);

    /* limiting case mu = Inf -> inverse chi-square as starting point*/
    if (!R_FINITE(mu))
	y = ACT_forceint(1/phi/qchisq(p, 1, !lower_tail, log_p));
    /* other cases -> Corning-Fisher */
    else
    {
	double sigma, sigma2, gamma;
	double phim2 = phi * mu * mu;

	sigma2 = phim2 * mu + mu;
	sigma = sqrt(sigma2);
	gamma = (3 * phim2 * sigma2 + mu)/sigma2/sigma;

	y = ACT_forceint(mu + sigma * (z + gamma * (z*z - 1)/6));
    }

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
