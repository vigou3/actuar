### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse pareto functions. The inverse Pareto
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = (x/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvpareto <- function (x, shape, scale, log = FALSE)
     .External("do_dpq", "dinvpareto", x, shape, scale, log)

pinvpareto <- function(q, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "pinvpareto", q, shape, scale, tail, log)

qinvpareto <- function(p, shape, scale, tail = TRUE, log = FALSE)
     .External("do_dpq", "qinvpareto", p, shape, scale, tail, log)

rinvpareto <- function(n, shape, scale)
     .External("do_random", "rinvpareto", n, shape, scale)
