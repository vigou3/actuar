### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}zmgeom functions to compute
### characteristics of the Zero Modified Geometric distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmgeom <- function (x, prob, p0, log = FALSE)
    .External(C_actuar_do_dpq, "dzmgeom", x, prob, p0, log)

pzmgeom <- function(q, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pzmgeom", q, prob, p0, lower.tail, log.p)

qzmgeom <- function(p, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qzmgeom", p, prob, p0, lower.tail, log.p)

rzmgeom <- function(n, prob, p0)
    .External(C_actuar_do_random, "rzmgeom", n, prob, p0)
