/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Function to compute the integral
 *
 *    G(a; x) = int_x^Infty t^(a - 1) exp(-t) dt
 *
 *  for a real and x > 0. When a > 0,
 *
 *    G(a; x) = gammafn(a) [1 - pgamma(x, a, scale = 1)].
 *
 *  When a < 0,
 *
 *    G(a; x) = -(x^a * exp(-x))/a + G(a + 1; x)/a,
 *            = (G(a + 1; x) - x^a * exp(-x))/a,
 *
 *  which can be repeated until a + k > 0. When a == 0,
 *
 *    G(0; x) = expint_E1(x),
 *
 *  the exponential integral.
 *
 *  See Appendix A of Klugman, Panjer & Willmot, Loss Models,
 *  Fourth Edition, Wiley, 2012 for the formula.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "dpq.h"
#include "actuar.h"

/* The workhorse */
double gammaint_raw(double x, double a)
{
/* Here, assume that (x, alph) are not NA and x > 0. */

     if (a > 0)
	return gammafn(a) *  pgamma(x, a, 1.0, /*l._t.*/0, /*log_p*/0);

    double k = ceil(-a);	/* note that k == 0 if a == 0 */
    double ap = a + k;		/* value a + k >= 0 */
    double lx = log(x);

    /* first computable value G(a + k; x) */
    double res = (ap == 0) ? expint_E1(x)
	: gammafn(ap) * pgamma(x, ap, 1.0, /*l._t.*/0, /*log_p*/0);

    /* Note that ap == 0 if a == 0 or if a == -1, -2, ... Enter the
     * recursive scheme only in the latter case, since we're done
     * otherwise. */
    while (ap > a)
    {
	ap--;
	res = (res - exp(ap * lx - x))/ap;
    }

    return res;
}

/* The function called by actuar_do_dpq(). The third argument is a
 * placeholder to fit into the scheme of dpq.c */
double gammaint(double x, double a, int foo)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a))
	return x + a;
#endif
    if (x <= 0)
    {
	if (x == 0 && a > 0)
	    return 1.0;
	else
	    return R_NaN;
    }

    return gammaint_raw(x, a);
}
