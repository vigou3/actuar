### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}zmpois functions to compute
### characteristics of the Zero Modified Poisson distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmpois <- function (x, lambda, p0, log = FALSE)
    .External(C_actuar_do_dpq, "dzmpois", x, lambda, p0, log)

pzmpois <- function(q, lambda, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pzmpois", q, lambda, p0, lower.tail, log.p)

qzmpois <- function(p, lambda, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qzmpois", p, lambda, p0, lower.tail, log.p)

rzmpois <- function(n, lambda, p0)
    .External(C_actuar_do_random, "rzmpois", n, lambda, p0)
