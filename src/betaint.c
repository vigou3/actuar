/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Function to compute the integral
 *
 *    B(a, b; x) = gammafn(a + b) int_0^x t^(a-1) (1-t)^(b-1) dt
 *
 *  for a > 0, b != -1, -2, ... and 0 < x < 1. When b > 0,
 *
 *    B(a, b; x) = gammafn(a) gammafn(b) pbeta(x, a, b).
 *
 *  When b < 0 and b != -1, -2, ... and a > 1 + floor(-b),
 *
 *    B(a, b; x)
 *    = -gammafn(a + b) {(x^(a-1) (1-x)^b)/b
 *        + [(a-1) x^(a-2) (1-x)^(b+1)]/[b(b+1)]
 *        + ...
 *        + [(a-1)...(a-r) x^(a-r-1) (1-x)^(b+r)]/[b(b+1)...(b+r)]}
 *      + [(a-1)...(a-r-1)]/[b(b+1)...(b+r)] gammafn(a-r-1)
 *        * gammafn(b+r+1) pbeta(x, a-r-1, b+r+1)
 *
 *  See Appendix A of Klugman, Panjer & Willmot, Loss Models,
 *  Fourth Edition, Wiley, 2012 for the formula.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "dpq.h"

double betaint_raw(double x, double a, double b)
{
/* Here, assume that (x, a, b) are not NA, 0 < x < 1 and a > 0. */

    if (b > 0)
	return gammafn(a) * gammafn(b) * pbeta(x, a, b, /*l._t.*/1, /*give_log*/0);

    double r = floor(-b);

    if (! (ACT_nonint(b) && a - r - 1 > 0))
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
     *  [previous term] * (a - i - 1)(1 - x)/[(b + i + 1) x]
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


/* The function called by actuar_do_dpq(). The fourth argument is a
 * placeholder to fit into the scheme of dpq.c */
double betaint(double x, double a, double b, int foo)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b))
	return x + a + b;
#endif
    if (a < 0 || x <= 0 || x >= 1)
	return R_NaN;

    return betaint_raw(x, a, b);
}



