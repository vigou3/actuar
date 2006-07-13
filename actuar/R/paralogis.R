### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}paralogistic functions. The paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^shape))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dparalogis <- function (x, shape, scale, log = FALSE)
  .External("do_dpq", "dparalogis", x, shape, scale, log)

pparalogis <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
  .External("do_dpq", "pparalogis", q, shape, scale, lower.tail, log.p)

qparalogis <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
  .External("do_dpq", "qparalogis", p, shape, scale, lower.tail, log.p)

rparalogis <- function(n, shape, scale)
  .External("do_random", "rparalogis", n, shape, scale)

mparalogis <- function(k, shape, scale, log = FALSE)
  .External("do_dpq", "mparalogis", k, shape, scale, log)

levparalogis <- function(x, shape, scale, order = 1, log = FALSE)
  .External("do_dpq", "levparalogis", x, shape, scale, order, log)
