### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}poisinvgauss functions to compute
### characteristics of the Poisson-Inverse Gaussian discrete
### distribution.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpoisinvgauss <- function(x, mean, shape = 1, dispersion = 1/shape, log = FALSE)
    .External(C_actuar_do_dpq, "dpoisinvgauss", x, mean, dispersion, log)

ppoisinvgauss <- function(q, mean, shape = 1, dispersion = 1/shape,
                          lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "ppoisinvgauss", q, mean, dispersion,
              lower.tail, log.p)

qpoisinvgauss <- function(p, mean, shape = 1, dispersion = 1/shape,
                          lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qpoisinvgauss", p, mean, dispersion,
              lower.tail, log.p)

rpoisinvgauss <- function(n, mean, shape = 1, dispersion = 1/shape)
    .External(C_actuar_do_random, "rpoisinvgauss", n, mean, dispersion)

## Aliases
dpig <- dpoisinvgauss
ppig <- ppoisinvgauss
qpig <- qpoisinvgauss
rpig <- rpoisinvgauss
