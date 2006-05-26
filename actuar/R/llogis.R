### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}loglogistic functions. The loglogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^shape/(1 + (x/scale)^shape)), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dllogis <- function (x, shape, scale, log = FALSE)
     .External("do_dpq", "dllogis", x, shape, scale, log)

pllogis <- function(q, shape, scale, tail = TRUE, log = FALSE)
     .External("do_dpq", "pllogis", q, shape, scale, tail, log)

qllogis <- function(p, shape, scale, tail = TRUE, log = FALSE)
     .External("do_dpq", "qllogis", p, shape, scale, tail, log)

rllogis <- function(n, shape, scale)
     .External("do_random", "rllogis", n, shape, scale)
