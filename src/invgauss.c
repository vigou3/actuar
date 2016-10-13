/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to calculate raw and limited moments for the inverse gaussian
 *  distribution. See ../R/InvGaussSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double dinvgauss(double x, double mu, double phi, int give_log)
{
    /*  We work with the density expressed as
     *
     *  (2 pi phi x^3)^(-1/2) exp(- u^2/(2 phi x))
     *
     *  with u = (x - mu)/mu.
     */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(phi))
	return x + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
	return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0)
	return R_FINITE(phi) ? ACT_D__0 : R_PosInf;

    /* limiting case phi = Inf and x > 0 */
    if (!R_FINITE(phi))
	return ACT_D__0;

    /* limiting case mu = Inf */
    if (!R_FINITE(mu))
	return ACT_D_exp(-(log(phi) + 3 * log(x) + 1/phi/x)/2 - M_LN_SQRT_2PI);

    /* standard cases */
    double xm = x/mu;
    double phim = phi * mu;

    return ACT_D_exp(-(log(phim) + 3 * log(xm) + R_pow_di(xm - 1, 2)/phim/xm)/2
		     - M_LN_SQRT_2PI - log(mu));
}

double pinvgauss(double q, double mu, double phi, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(q) || ISNAN(mu) || ISNAN(phi))
	return q + mu + phi;
#endif
    if (mu <= 0.0 || phi <= 0.0)
        return R_NaN;

    if (q < 0)
        return ACT_DT_0;

    /* handle x == 0 separately */
    if (q == 0)
	return R_FINITE(phi) ? ACT_DT_0 : ACT_DT_1;

    /* limiting cases phi = Inf and q > 0, and q = Inf */
    if (!R_FINITE(phi) || !R_FINITE(q))
	return ACT_DT_1;

    /* limiting case mu = Inf */
    if (!R_FINITE(mu))
	return pchisq(1/q/phi, 1, !lower_tail, log_p);

    /* standard cases */
    double qm = q/mu;
    double phim = phi * mu;

    /* approximation for (survival) probabilities in the far right tail */
    if (!lower_tail && qm > 1e6)
    {
	double r = qm/2/phim;
	if (r > 5e5)
	    return ACT_D_exp(1/phim - M_LN_SQRT_PI - log(2*phim) - 1.5 * log1p(r) - r);
    }

    /* all other probabilities rely on pnorm() */
    double a, b, r = sqrt(q * phi);

    a = pnorm((qm - 1)/r, 0, 1, lower_tail, /* log_p */1);
    b = 2/phim + pnorm(-(qm + 1)/r, 0, 1, /* l._t. */1, /* log_p */1);

    return ACT_D_exp(a + (lower_tail ? log1p(exp(b - a)) : ACT_Log1_Exp(b - a)));
}


double minvGauss(double order, double nu, double lambda, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(nu) || ISNAN(lambda))
	return order + nu + lambda;
#endif
    if (!R_FINITE(nu) ||
        !R_FINITE(lambda) ||
        !R_FINITE(order) ||
        nu <= 0.0 ||
        lambda <= 0.0 ||
        (int) order != order)
        return R_NaN;

    /* Trivial case */
    if (order == 0.0)
        return 0.0;

    int i, n = order;
    double z = 0.0;

    for (i = 0; i < n; i++)
        z += R_pow_di(nu, n) * gammafn(n + i) *
            R_pow_di(2.0 * lambda/nu, -i) /
            (gammafn(i + 1) * gammafn(n - i));
    return z;
}

double levinvGauss(double limit, double nu, double lambda, double order,
                   int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(nu) || ISNAN(lambda) || ISNAN(order))
	return limit + nu + lambda + order;
#endif
    if (!R_FINITE(nu)     ||
        !R_FINITE(lambda) ||
        !R_FINITE(order)  ||
        nu <= 0.0    ||
        lambda < 0.0 ||
        order != 1.0)
        return R_NaN;

    if (limit <= 0.0)
        return 0.0;

    /* From R, order == 1 */
    double tmp, y, z;

    tmp = sqrt(lambda/limit);
    y = (limit + nu)/nu;
    z = (limit - nu)/nu;

    return limit - nu * z * pnorm(z * tmp, 0.0, 1.0, 1, 0)
        - nu * y * exp(2.0 * lambda/nu) * pnorm(-y * tmp, 0.0, 1.0, 1, 0);
}

double mgfinvGauss(double x, double nu, double lambda, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(nu) || ISNAN(lambda))
	return x + nu + lambda;
#endif
    if (!R_FINITE(nu) ||
        !R_FINITE(lambda) ||
        nu <= 0.0 ||
        lambda < 0.0 ||
        x > lambda/(2.0 * nu * nu))
        return R_NaN;

    if (x == 0.0)
        return ACT_D_exp(0.0);

    return ACT_D_exp(lambda / nu * (1.0 - sqrt(1.0 - 2.0 * nu * nu * x/lambda)));
}
