### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev,mgf}invgauss functions to compute
### raw and limited moments, and the moment generating function for
### the Inverse Gaussian distribution.
###
### See Part 2 of Chhikara and Folks, The inverse Gaussian
### distribution: theory, methodology and application, Decker, 1989;
### Part 1 of Seshadri, The inverse Gaussian distribution: statistical
### theory and applications, Springer, 1989
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvgauss <- function(x, mean, shape = 1, dispersion = 1/shape, log = FALSE)
    .External("actuar_do_dpq", "dinvgauss", x, mean, dispersion, log)

pinvgauss <- function(q, mean, shape = 1, dispersion = 1/shape,
                      lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvgauss", q, mean, dispersion,
              lower.tail, log.p)

qinvgauss <- function(p, mean, shape = 1, dispersion = 1/shape,
                      lower.tail = TRUE, log.p = FALSE,
                      tol = 1e-14, maxit = 100, echo = FALSE, trace = echo)
    .External("actuar_do_dpq", "qinvgauss", p, mean, dispersion,
              lower.tail, log.p, tol, maxit, trace)

rinvgauss <- function(n, mean, shape = 1, dispersion = 1/shape)
    .External("actuar_do_random", "rinvgauss", n, mean, dispersion)

minvgauss <- function(order, mean, shape = 1, dispersion = 1/shape)
    .External("actuar_do_dpq", "minvgauss", order, mean, dispersion, FALSE)

levinvGauss <- function(limit, nu, lambda, order = 1)
    .External("actuar_do_dpq", "levinvGauss", limit,  nu, lambda, order, FALSE)

mgfinvGauss <- function(x, nu, lambda, log = FALSE)
    .External("actuar_do_dpq", "mgfinvGauss", x, nu, lambda, log)

## Functions deprecated in actuar v2.0-0
minvGauss <- function(order, nu, lambda)
{
    .Deprecated("minvgauss", package = "actuar")
    .External("actuar_do_dpq", "minvGauss", order, nu, lambda, FALSE)
}

## aliases
levinvgauss <- levinvGauss
mgfinvgauss <- mgfinvGauss
