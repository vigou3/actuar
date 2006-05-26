### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse paralogistic functions. The inverse paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^tau/(1 + (x/scale)^tau))^tau, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvparalogis <- function (x, shape, scale, log = FALSE)
    .External("do_dpq", "dinvparalogis", x, shape, scale, log)

pinvparalogis <- function(q, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "pinvparalogis", q, shape, scale, tail, log)

qinvparalogis <- function(p, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "qinvparalogis", p, shape, scale, tail, log)

rinvparalogis <- function(n, shape, scale)
    .External("do_random", "rinvparalogis", n, shape, scale)
