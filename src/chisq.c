/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Functions to calculate raw and limited moments for the Chi-square
 *  distribution. See ../R/ChisqSupp.R for details.
 *
 *  AUTHORS: Christophe Dutang and Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "locale.h"
#include "dpq.h"

double mchisq(double order, double df, double ncp, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(order) || ISNAN(df) || ISNAN(ncp))
	return order + df + ncp;
#endif
    if (!R_FINITE(df)    ||
        !R_FINITE(ncp)   ||
        !R_FINITE(order) ||
        df <= 0.0 ||
        ncp < 0.0)
        return R_NaN;

    if (order <= -df/2.0)
	return R_PosInf;

    /* Trivial case */
    if (order == 0.0)
        return 1.0;

    /* Centered chi-square distribution */
    if (ncp == 0.0)
        return R_pow(2.0, order) * gammafn(order + df/2.0) / gammafn(df/2.0);

    /* Non centered chi-square distribution */
    if (order >= 1.0 && (int) order == order)
    {
        int i, j = 0, n = order;
        double *res;

        /* Array with 1, E[X], E[X^2], ..., E[X^n] */
        res = (double *) R_alloc(n + 1, sizeof(double));
        res[0] = 1.0;
        res[1] = df + ncp;      /* E[X] */
        for (i = 2; i <= n; i++)
        {
            res[i] = R_pow_di(2.0, i - 1) * (df + i * ncp);
            for (j = 1; j < i; j++)
                res[i] += R_pow_di(2.0, j - 1) * (df + j * ncp) * res[i - j] / gammafn(i - j + 1);
            res[i] *= gammafn(i);
        }
        return res[n];
    }
    else
        return R_NaN;
}

double levchisq(double limit, double df, double ncp, double order, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(limit) || ISNAN(df) || ISNAN(ncp) || ISNAN(order))
	return limit + df + ncp + order;
#endif
    if (!R_FINITE(df)    ||
        !R_FINITE(ncp)   ||
        !R_FINITE(order) ||
        df <= 0.0 ||
        ncp < 0.0)
        return R_NaN;

    if (order <= -df/2.0)
	return R_PosInf;

    if (limit <= 0.0)
        return 0.0;

    if (ncp == 0.0)
    {
        double u, tmp;

        tmp = order + df/2.0;
        u = exp(log(limit) - M_LN2);

        return R_pow(2.0, order) * gammafn(tmp) *
            pgamma(u, tmp, 1.0, 1, 0) / gammafn(df/2.0) +
            ACT_DLIM__0(limit, order) * pgamma(u, df/2.0, 1.0, 0, 0);
    }
    else
        return R_NaN;
}

double mgfchisq(double t, double df, double ncp, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(t) || ISNAN(df) || ISNAN(ncp))
	return t + df + ncp;
#endif
    if (!R_FINITE(df)  ||
        !R_FINITE(ncp) ||
        df <= 0.0 ||
        ncp < 0.0 ||
        2.0 * t > 1.0)
        return R_NaN;

    if (t == 0.0)
        return ACT_D__1;

    return ACT_D_exp(ncp * t / (1.0 - 2.0 * t) - df/2.0 * log1p(-2.0 * t));
}
