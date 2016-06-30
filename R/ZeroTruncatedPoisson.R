### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r}ztpois functions to compute
### characteristics of the Zero Truncated Poisson distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dztpois <- function (x, lambda, log = FALSE)
    .External("actuar_do_dpq", "dztpois", x, lambda, log)

pztpois <- function(q, lambda, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pztpois", q, lambda, lower.tail, log.p)

qztpois <- function(p, lambda, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qztpois", p, lambda, lower.tail, log.p)

rztpois <- function(n, lambda)
    .External("actuar_do_random", "rztpois", n, lambda)
