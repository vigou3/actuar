### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}paralogistic functions. The paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^alpha))^alpha, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dparalogis <- function (x, shape, scale, log = FALSE)
    .External("do_dpq", "dparalogis", x, shape, scale, log)

#pparalogis <- function(q, shape, scale, tail = TRUE, log = FALSE)
#     .External("do_dpq", "pparalogis", q, shape, scale, tail, log)

#qparalogis <- function(p, shape, scale, tail = TRUE, log = FALSE)
#     .External("do_dpq", "qparalogis", p, shape, scale, tail, log)

rparalogis <- function(n, shape, scale)
    .External("do_random", "rparalogis", n, shape, scale)
