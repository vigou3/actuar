### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse paralogistic functions. The inverse paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^shape/(1 + (x/scale)^shape))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvparalogis <- function (x, shape, scale, log = FALSE)
    .External("do_dpq", "dinvparalogis", x, shape, scale, log)

pinvparalogis <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "pinvparalogis", q, shape, scale, lower.tail, log.p)

qinvparalogis <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qinvparalogis", p, shape, scale, lower.tail, log.p)

rinvparalogis <- function(n, shape, scale)
    .External("do_random", "rinvparalogis", n, shape, scale)
