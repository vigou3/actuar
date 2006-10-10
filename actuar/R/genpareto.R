### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}genpareto functions to compute
### characteristics of the Generalized Pareto distribution.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dgenpareto <- function(x, shape1, shape2, rate = 1, scale = 1/rate,
                       log = FALSE)
     .External("do_dpq", "dgenpareto", x, shape1, shape2, scale, log)

pgenpareto <- function(q, shape1, shape2, rate = 1, scale = 1/rate,
                       lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "pgenpareto", q, shape1, shape2, scale,
               lower.tail, log.p)

qgenpareto <- function(p, shape1, shape2, rate = 1, scale = 1/rate,
                       lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qgenpareto", p, shape1, shape2, scale,
               lower.tail, log.p)

rgenpareto <- function(n, shape1, shape2, rate = 1, scale = 1/rate)
     .External("do_random", "rgenpareto", n, shape1, shape2, scale)

mgenpareto <- function(order, shape1, shape2, rate = 1, scale = 1/rate,
                       log = FALSE)
     .External("do_dpq", "mgenpareto", order, shape1, shape2, scale, log)

levgenpareto <- function(limit, shape1, shape2, rate = 1, scale = 1/rate,
                         order = 1, log = FALSE)
     .External("do_dpq", "levgenpareto", limit, shape1, shape2, scale,
               order, log)
