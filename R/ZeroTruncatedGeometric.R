### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}ztgeom functions to compute
### characteristics of the Zero Truncated Geometric distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dztgeom <- function(x, prob, log = FALSE)
    .External(C_actuar_do_dpq, "dztgeom", x, prob, log)

pztgeom <- function(q, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pztgeom", q, prob, lower.tail, log.p)

qztgeom <- function(p, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qztgeom", p, prob, lower.tail, log.p)

rztgeom <- function(n, prob)
    .External(C_actuar_do_random, "rztgeom", n, prob)
