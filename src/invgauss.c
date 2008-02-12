/*  ===== actuar: an R package for Actuarial Science =====
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

double minvGauss(double order, double nu, double lambda, int give_log)
{
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
            R_pow_di(2 * lambda/nu, -i) /
            (gammafn(i + 1) * gammafn(n - i));
    return z;
}

double levinvGauss(double limit, double nu, double lambda, double order,
                   int give_log)
{
    double tmp, y, z;

    if (!R_FINITE(nu)     ||
        !R_FINITE(lambda) ||
        !R_FINITE(order)  ||
        nu <= 0.0    ||
        lambda < 0.0 ||
        order != 1.0)
        return R_NaN;

    if (limit <= 0.0)
        return 0;

    /* From R, order == 1 */
    tmp = sqrt(lambda/limit);
    y = (limit + nu)/nu;
    z = (limit - nu)/nu;

    return limit - nu * z * pnorm(z * tmp, 0.0, 1.0, 0, 0)
        - nu * y * exp(2 * lambda/nu) * pnorm(-y * tmp, 0., 1.0, 0, 0);
}

double mgfinvGauss(double x, double nu, double lambda, int give_log)
{
    if (!R_FINITE(nu) ||
        !R_FINITE(lambda) ||
        nu <= 0.0 ||
        lambda < 0.0 ||
        x > lambda/(2 * nu * nu))
        return R_NaN;

    if (x == 0.0)
        return R_D_exp(0.0);

    return R_D_exp(lambda / nu * (1 - sqrt(1 - 2 * nu * nu * x/lambda)));
}
