### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r}zmnbinom functions to compute
### characteristics of the Zero Modified Negative Binomial
### distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmnbinom <- function (x, size, prob, p0, log = FALSE)
    .External("actuar_do_dpq", "dzmnbinom", x, size, prob, p0, log)

pzmnbinom <- function(q, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pzmnbinom", q, size, prob, p0, lower.tail, log.p)

qzmnbinom <- function(p, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qzmnbinom", p, size, prob, p0, lower.tail, log.p)

rzmnbinom <- function(n, size, prob, p0)
    .External("actuar_do_random", "rzmnbinom", n, size, prob, p0)

rzmnbinom2 <- function(n, size, prob, p0)
    .External("actuar_do_random", "rzmnbinom2", n, size, prob, p0)
