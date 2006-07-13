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

pllogis <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "pllogis", q, shape, scale, lower.tail, log.p)

qllogis <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qllogis", p, shape, scale, lower.tail, log.p)

rllogis <- function(n, shape, scale)
     .External("do_random", "rllogis", n, shape, scale)

mllogis <- function(k, shape, scale, log = FALSE)
     .External("do_dpq", "mllogis", k, shape, scale, log)

levllogis <- function(x, shape, scale, order = 1, log = FALSE)
     .External("do_dpq", "levllogis", x, shape, scale, order, log)
