/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and quantile
 *  functions, raw and limited moments and to simulate random variates
 *  for the transformed beta distribution. See ../R/TransformedBeta.R for
 *  details.
 *
 *  AUTHORS: Mathieu Pigeon and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

/* Function to compute a scaled version of the beta distribution
 * function when the second parameter is a negative non integer. More
 * precisely, the function computes the value of
 *
 *   gamma(a) gamma(b) pbeta(x, a, b)
 *
 * when b < 0, b is not an integer, and a > 1 + floor(-b).
 *
 * See Appendix A of Klugman, Panjer & Willmot, Loss Models,
 * Second Edition, Wiley, 2004 for the formula.
 */

double actuar_pbetanegb(double x, double a, double b)
{
    double r = floor(-b);

    /* Return NaN if b is integer or if a is too small for the value
     * of b. */
    if (b == (int) b || a - r - 1 <= 0)
	return R_NaN;

    /* Most of the effort is spent on computing the sum that is
     * alternating since b < 0 in the denominators. The terms are
     * calculated in the log scale. Furthermore, negative and
     * positive terms are accumulated in separate variables and the
     * sum (the difference, actually) is done at the end. */
    int i, s;
    double ap = a, bp = b;		/* copies of a and b */
    double res1, res2, lasum, lbsum;
    double lx = log(x), lx1m = log1p(-x);

    /* lasum and lbsum will hold log[(a - 1)(a - 2) ... (a - r)] and
     * log[(-b)(-b - 1)...(-b - r)], respectively */
    lasum = 0.;			/* initialize to 0 */
    lbsum = log(-bp);		/* initialize to log(-b) */

    /* Computation of the first term in the alternating sum. Sign is
     * negative. */
    ap--;			/* a - 1 */
    res1 = exp(ap * lx + bp * lx1m - lbsum); /* holds the sum of negative terms */
    bp++;			/* b + 1 */

    /* Simplest case with only the first term. */
    if (r == 0.)
 	return(-gammafn(a + b) * (-res1)
	       - exp(log(ap) + lgammafn(ap) - lbsum + lgammafn(bp) + pbeta(x, ap, bp, 1, 1)));

    /* Other terms in the alternating sum iff r > 0. */
    res2 = 0.;			/* holds the sum of positive terms */
    s = -1;			/* sign of computed term; previous one was negative */

    for (i = 0; i < r; i++)
    {
	s = -s;			/* sign change */
	lasum += log(ap);	/* log(a - 1) + ... + log(a - i - 1) */
	lbsum += log(-bp);	/* log(-b) + ... + log(-b - i - 1) */
	ap--;
	if (s > 0)		/* positive terms */
	    res2 += exp(lasum + ap * lx + bp * lx1m - lbsum);
	else			/* negative terms */
	    res1 += exp(lasum + ap * lx + bp * lx1m - lbsum);
	bp++;
    }

    return(-gammafn(a + b) * (res2 - res1) +
	   (double) s * exp(lasum + log(ap) + lgammafn(ap) - lbsum + lgammafn(bp) + pbeta(x, ap, bp, 1, 1)));
}

double actuar_pbetanegb1(double x, double a, double b)
{
    double r = floor(-b);

    /* Return NaN if b is integer or if a is too small for the value
     * of b. */
    if (b == (int) b || a - r - 1 <= 0)
	return R_NaN;

    /* Most of the effort is spent on computing the sum that is
     * alternating since b < 0 in the denominators. The terms are
     * calculated in the log scale, except for the denominators [b(b +
     * 1) ... (b + r)]. This product is computed (and accumulated)
     * separately to 1) avoid changing the sign of b to work with
     * logs; 2) take care of the sign of the product automatically. */
    int i;
    double ap = a, bp = b;		/* copies of a and b */
    double sum, lasum, bprod;
    double lx = log(x), lx1m = log1p(-x);

    /* Computation of the first term in the alternating sum. */
    ap--;				 /* a - 1 */
    sum = exp(ap * lx + bp * lx1m) / bp; /* holds the sum of negative terms */
    bp++;				 /* b + 1 */

    lasum = 0.;
    bprod = bp;

    /* Other terms in the alternating sum iff r > 0. */
    for (i = 0; i < r; i++)
    {
	lasum += log(ap);	/* log(a - 1) + ... + log(a - i - 1) */
	bprod *= bp;		/* b (b + 1) ... (b + i + 1) */
	ap--;
	sum += exp(lasum + ap * lx + bp * lx1m) / bprod;
	bp++;
    }

    return(-gammafn(a + b) * sum +
	   exp(lasum + log(ap) + lgammafn(ap) + lgammafn(bp) + pbeta(x, ap, bp, 1, 1)) / bprod);
}


double actuar_pbetanegb2(double x, double a, double b)
{
    double r = floor(-b);

    /* Return NaN if b is integer or if a is too small for the value
     * of b. */
    if (b == (int) b || a - r - 1 <= 0)
	return R_NaN;

    /* There are two quantities to accumulate in order to compute the
     * final result: the alternating sum (to be stored in 'sum') and
     * the ratio [(a - 1) ... (a - r)]/[b(b + 1) ... (b + r)] (to be
     * stored in 'ratio'). Some calculations are done in the log
     * scale. */
    int i;
    double ap = a, bp = b;		/* copies of a and b */
    double c, sum, ratio;
    double lx = log(x), lx1m = log1p(-x), x1 = exp(lx1m - lx);

    /* Computation of the first term in the alternating sum. */
    ap--;			       /* a - 1 */
    c = exp(ap * lx + bp * lx1m) / bp; /* (x^(a - 1) (1 - x)^b) / b */
    sum = c;			       /* first term */
    ratio = 1. / bp;		       /* 1 / b */
    bp++;			       /* b + 1 */

    /* Other terms in the alternating sum iff r > 0. */
    for (i = 0; i < r; i++)
    {
	c *= (ap / bp) * x1;	       /* new term in the sum */
	sum += c;
	ratio *= ap / bp;
	ap--;
	bp++;
    }

    return(-gammafn(a + b) * sum
	   + (ratio * ap) * exp(lgammafn(ap) + lgammafn(bp) + pbeta(x, ap, bp, 1, 1)));
}


double dtrbeta(double x, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
    /*  We work with the density expressed as
     *
     *  shape2 * u^shape3 * (1 - u)^shape1 / (x * beta(shape1, shape3))
     *
     *  with u = v/(1 + v) = 1/(1 + 1/v), v = (x/scale)^shape2.
     */

    double tmp, logu, log1mu;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    if (!R_FINITE(x) || x < 0.0)
        return ACT_D__0;

    /* handle x == 0 separately */
    if (x == 0.0)
    {
	if (shape2 * shape3 < 1) return R_PosInf;
	if (shape2 * shape3 > 1) return ACT_D__0;
	/* else */
	return give_log ?
	    log(shape2) - log(scale) - lbeta(shape3, shape1) :
	    shape2 / (scale * beta(shape3, shape1));
    }

    tmp = shape2 * (log(x) - log(scale));
    logu = - log1pexp(-tmp);
    log1mu = - log1pexp(tmp);

    return ACT_D_exp(log(shape2) + shape3 * logu + shape1 * log1mu
                   - log(x) - lbeta(shape3, shape1));
}

double ptrbeta(double q, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
    double u;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (q <= 0)
        return ACT_DT_0;

    u = exp(-log1pexp(-shape2 * (log(q) - log(scale))));

    return pbeta(u, shape3, shape1, lower_tail, log_p);
}

double qtrbeta(double p, double shape1, double shape2, double shape3,
               double scale, int lower_tail, int log_p)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    ACT_Q_P01_boundaries(p, 0, R_PosInf);
    p = ACT_D_qIv(p);

    return scale * R_pow(1.0 / qbeta(p, shape3, shape1, lower_tail, 0) - 1.0,
                         -1.0 / shape2);
}

double rtrbeta(double shape1, double shape2, double shape3, double scale)
{
    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale <= 0.0)
        return R_NaN;

    return scale * R_pow(1.0 / rbeta(shape3, shape1) - 1.0, -1.0 / shape2);
}

double mtrbeta(double order, double shape1, double shape2, double shape3,
               double scale, int give_log)
{
    double tmp;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale)  ||
        !R_FINITE(order)  ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
	return R_NaN;

    if (order <= - shape3 * shape2 ||
        order >= shape1 * shape2)
        return R_PosInf;

    tmp = order / shape2;

    return R_pow(scale, order) * beta(shape3 + tmp, shape1 - tmp)
        / beta(shape1, shape3);
}

double levtrbeta(double limit, double shape1, double shape2, double shape3,
                 double scale, double order, int give_log)
{
    double u, tmp, a, b, r;

    if (!R_FINITE(shape1) ||
        !R_FINITE(shape2) ||
        !R_FINITE(shape3) ||
        !R_FINITE(scale) ||
        !R_FINITE(order) ||
        shape1 <= 0.0 ||
        shape2 <= 0.0 ||
        shape3 <= 0.0 ||
        scale  <= 0.0)
        return R_NaN;

    if (order  <= - shape3 * shape2)
        return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    r = order / shape2;
    a = shape3 + r;
    b = shape1 - r;

    u = exp(-log1pexp(-shape2 * (log(limit) - log(scale))));

    if (b < 0)
	tmp = give_log ? actuar_pbetanegb2(u, a, b) / gammafn(shape1 + shape3)
	    : actuar_pbetanegb(u, a, b) / gammafn(shape1 + shape3);
    else
	tmp = beta(a, b) * pbeta(u, a, b, 1, 0);

    return R_pow(scale, order) * tmp / beta(shape1, shape3)
	+ ACT_DLIM__0(limit, order) * pbeta(u, shape3, shape1, 0, 0);
}
