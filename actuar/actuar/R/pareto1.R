### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}single-parameter pareto functions. The single-parameter ### Pareto distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (min/x)^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto1 <- function (x, shape, min, log = FALSE)
    .External("do_dpq", "dpareto1", x, shape, min, log)

ppareto1 <- function(q, shape, min, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "ppareto1", q, shape, min, lower.tail, log.p)

qpareto1 <- function(p, shape, min, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qpareto1", p, shape, min, lower.tail, log.p)

rpareto1 <- function(n, shape, min)
    .External("do_random", "rpareto1", n, shape, min)
