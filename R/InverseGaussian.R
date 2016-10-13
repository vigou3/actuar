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

