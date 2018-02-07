### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}zmbinom functions to compute
### characteristics of the Zero Modified Binomial distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmbinom <- function (x, size, prob, p0, log = FALSE)
    .External(C_actuar_do_dpq, "dzmbinom", x, size, prob, p0, log)

pzmbinom <- function(q, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pzmbinom", q, size, prob, p0, lower.tail, log.p)

qzmbinom <- function(p, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qzmbinom", p, size, prob, p0, lower.tail, log.p)

rzmbinom <- function(n, size, prob, p0)
    .External(C_actuar_do_random, "rzmbinom", n, size, prob, p0)
