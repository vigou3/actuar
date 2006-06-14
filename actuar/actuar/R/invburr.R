### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse burr functions. The inverse Burr
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^shape2/(1 + (x/scale)^shape2))^shape1 x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvburr <- function (x, shape1, scale, shape2, log = FALSE)
    .External("do_dpq", "dinvburr", x, shape1, scale, shape2, log)

pinvburr <- function(q, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "pinvburr", q, shape1, scale, shape2, lower.tail, log.p)

qinvburr <- function(p, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qinvburr", p, shape1, scale, shape2, lower.tail, log.p)

rinvburr <- function(n, shape1, scale, shape2)
    .External("do_random", "rinvburr", n, shape1, scale, shape2)
