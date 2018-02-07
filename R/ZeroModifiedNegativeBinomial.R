### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}zmnbinom functions to compute
### characteristics of the Zero Modified Negative Binomial
### distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmnbinom <- function (x, size, prob, p0, log = FALSE)
    .External(C_actuar_do_dpq, "dzmnbinom", x, size, prob, p0, log)

pzmnbinom <- function(q, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pzmnbinom", q, size, prob, p0, lower.tail, log.p)

qzmnbinom <- function(p, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qzmnbinom", p, size, prob, p0, lower.tail, log.p)

rzmnbinom <- function(n, size, prob, p0)
    .External(C_actuar_do_random, "rzmnbinom", n, size, prob, p0)
