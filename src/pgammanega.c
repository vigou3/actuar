/*  ===== actuar: An R Package for Actuarial Science =====
 *
 *  Function to compute the integral
 *
 *    G(a; x) = gammafn(a) [1 - pgamma(x, a, scale = 1)]
 *
 *  where a *may* be nonpositive. If a == 0,
 *
 *    G(0; x) = expint_E1(x).
 *
 *  If a < 0, we have
 *
 *    G(a; x) = - (x^a * exp(-x))/a + G(a + 1; x)/a,
 *
 *  which can be repeated until a + k > 0.
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

double pgammanega_scaled(double x, double a,
			 int foo /*unused; for scheme of dpq.c */)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a))
	return x + a;
#endif

    /* simple case using pgamma() directly */
    if (a > 0)
	return gammafn(a) *  pgamma(x, a, 1.0, /*l._t.*/0, /*log_p*/0);

    double k = ceil(-a);	/* note that k == 0 if a == 0 */
    double ap = a + k;		/* value a + k >= 0 */

    /* first computable value G(a + k; x) */
    double res = (ap == 0) ? expint_E1(x, 0.0, 0)
	: gammafn(ap) * pgamma(x, ap, 1.0, /*l._t.*/0, /*log_p*/0);

    double lx = log(x);

    /* note that ap == 0 if a is zero or negative integer */

    while (ap > a)
    {
	ap--;
	res = (res - exp(ap * lx - x))/ap;
    }

    Rprintf("G(a; x) = %f\n", res);
    return res;
}
