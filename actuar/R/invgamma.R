### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}inverse gamma functions. The inverse gamma
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - pgamma(scale/x, shape, scale), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvgamma <- function (x, shape, scale, log = FALSE)
     .External("do_dpq", "dinvgamma", x, shape, scale, log)

pinvgamma <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "pinvgamma", q, shape, scale, lower.tail, log.p)

qinvgamma <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qinvgamma", p, shape, scale, lower.tail, log.p)

rinvgamma <- function(n, shape, scale)
     .External("do_random", "rinvgamma", n, shape, scale)

minvgamma <- function(k, shape, scale, log = FALSE)
     .External("do_dpq", "minvgamma", k, shape, scale, log)

levinvgamma <- function(d, shape, scale, order = 1, log = FALSE)
     .External("do_dpq", "levinvgamma", d, shape, scale, order, log)
