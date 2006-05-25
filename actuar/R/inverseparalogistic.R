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

diparalogis <- function (x, shape, scale, log = FALSE)
    .External("do_dpq", "diparalogis", x, shape, scale, log)

##piparalogis <- function(q, shape, scale, tail = TRUE, log = FALSE)
##    .External("do_dpq", "piparalogis", q, shape, scale, tail, log)

##qiparalogis <- function(p, shape, scale, tail = TRUE, log = FALSE)
##    .External("do_dpq", "qiparalogis", p, shape, scale, tail, log)

riparalogis <- function(n, shape, scale)
    .External("do_random", "riparalogis", n, shape, scale)
