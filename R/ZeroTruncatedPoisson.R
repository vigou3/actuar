### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}ztpois functions to compute
### characteristics of the Zero Truncated Poisson distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dztpois <- function (x, lambda, log = FALSE)
    .External(C_actuar_do_dpq, "dztpois", x, lambda, log)

pztpois <- function(q, lambda, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pztpois", q, lambda, lower.tail, log.p)

qztpois <- function(p, lambda, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qztpois", p, lambda, lower.tail, log.p)

rztpois <- function(n, lambda)
    .External(C_actuar_do_random, "rztpois", n, lambda)

## not exported; for internal use in panjer()
pgfztpois <- function(x, lambda)
    expm1(lambda * x)/expm1(lambda)
