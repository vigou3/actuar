/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Function to compute a scaled version of the beta distribution
 *  function when the second parameter is a negative non integer. More
 *  precisely, the function computes the value of
 *
 *    gamma(a) gamma(b) pbeta(x, a, b)
 *
 *  when b < 0, b != -1, -2, ... and a > 1 + floor(-b).
 *
 *  See Appendix A of Klugman, Panjer & Willmot, Loss Models,
 *  Fourth Edition, Wiley, 2012 for the formula.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "dpq.h"

double pbetanegb(double x, double a, double b, int foo /*unused*/)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b))
	return x + a + b;
#endif

    double r = floor(-b);

    /* Return NaN if b is integer or if a is too small for the value
     * of b. */
    if (!(ACT_nonint(b) && a - r - 1 > 0))
	return R_NaN;

    /* There are two quantities to accumulate in order to compute the
     * final result: the alternating sum (to be stored in 'sum') and
     * the ratio [(a - 1) ... (a - r)]/[b(b + 1) ... (b + r)] (to be
     * stored in 'ratio'). Some calculations are done in the log
     * scale. */
    int i;
    double ap = a, bp = b;		/* copies of a and b */
    double lx = log(x);	                /* log(x) */
    double lx1m = log1p(-x);		/* log(1 - x) */;
    double x1 = exp(lx1m - lx);         /* (1 - x)/x */
    double c, tmp, sum, ratio;

    /* Computation of the first term in the alternating sum. */
    ap--;			       /* a - 1 */
    c = exp(ap * lx + bp * lx1m)/bp;   /* (x^(a - 1) (1 - x)^b) / b */
    sum = c;			       /* first term */
    ratio = 1/bp;		       /* 1 / b */
    bp++;			       /* b + 1 */

    /* Other terms in the alternating sum iff r > 0.
     * Relies on the fact that each new term in the sum is
     *
     *  previous term * (a - i - 1)(1 - x)/[(b + i + 1) x]
     *
     * for i = 0, ..., r - 1
     */
    for (i = 0; i < r; i++)
    {
	tmp = ap/bp;		/* (a - i - 1)/(b + i + 1) */
	c *= tmp * x1;		/* new term in the sum  */
	sum += c;
	ratio *= tmp;
	ap--;
	bp++;
    }

    return(-gammafn(a + b) * sum
	   + (ratio * ap) * exp(lgammafn(ap) + lgammafn(bp) +
				pbeta(x, ap, bp, /*l._t.*/1, /*give_log*/1)));
}
