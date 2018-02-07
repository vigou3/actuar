### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev,mgf}invgauss functions to compute
### characteristics of the Inverse Gaussian distribution.
###
### Functions [dpq]invgauss rely on C implementations of functions of
### the same name in package statmod. See:
###
### Giner, G. and Smyth, G. K. (2016), "statmod: Probability
###   Calculations for the Inverse Gaussian Distribution", R Journal,
###   vol. 8, no 1, p. 339-351.
###   https://journal.r-project.org/archive/2016-1/giner-smyth.pdf
###
### Chhikara, R. S. and Folk, T. L. (1989), The Inverse Gaussian
###   Distribution: Theory, Methodology and Applications}, Decker.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvgauss <- function(x, mean, shape = 1, dispersion = 1/shape, log = FALSE)
    .External(C_actuar_do_dpq, "dinvgauss", x, mean, dispersion, log)

pinvgauss <- function(q, mean, shape = 1, dispersion = 1/shape,
                      lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pinvgauss", q, mean, dispersion,
              lower.tail, log.p)

qinvgauss <- function(p, mean, shape = 1, dispersion = 1/shape,
                      lower.tail = TRUE, log.p = FALSE,
                      tol = 1e-14, maxit = 100, echo = FALSE, trace = echo)
    .External(C_actuar_do_dpq, "qinvgauss", p, mean, dispersion,
              lower.tail, log.p, tol, maxit, trace)

rinvgauss <- function(n, mean, shape = 1, dispersion = 1/shape)
    .External(C_actuar_do_random, "rinvgauss", n, mean, dispersion)

minvgauss <- function(order, mean, shape = 1, dispersion = 1/shape)
    .External(C_actuar_do_dpq, "minvgauss", order, mean, dispersion, FALSE)

levinvgauss <- function(limit, mean, shape = 1, dispersion = 1/shape, order = 1)
    .External(C_actuar_do_dpq, "levinvgauss", limit, mean, dispersion, order, FALSE)

mgfinvgauss <- function(t, mean, shape = 1, dispersion = 1/shape, log = FALSE)
    .External(C_actuar_do_dpq, "mgfinvgauss", t, mean, dispersion, log)

## Functions deprecated in actuar v2.0-0
minvGauss <- function(order, nu, lambda)
{
    .Deprecated("minvgauss", package = "actuar")
    .External(C_actuar_do_dpq, "minvGauss", order, nu, lambda, FALSE)
}

levinvGauss <- function(limit, nu, lambda, order = 1)
{
    .Deprecated("levinvgauss", package = "actuar")
    .External(C_actuar_do_dpq, "levinvGauss", limit,  nu, lambda, order, FALSE)
}

mgfinvGauss <- function(x, nu, lambda, log = FALSE)
{
    .Deprecated("mgfinvgauss", package = "actuar")
    .External(C_actuar_do_dpq, "mgfinvGauss", x, nu, lambda, log)
}
