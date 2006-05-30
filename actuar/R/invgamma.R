### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse gamma functions. The inverse gamma
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

pinvgamma <- function(q, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "pinvgamma", q, shape, scale, tail, log)

qinvgamma <- function(p, shape, scale, tail = TRUE, log = FALSE)
     .External("do_dpq", "qinvgamma", p, shape, scale, tail, log)

rinvgamma <- function(n, shape, scale)
     .External("do_random", "rinvgamma", n, shape, scale)
