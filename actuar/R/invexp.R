### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse exponential functions. The
### inverse exponential distribution used in these functions has
### cumulative distribution function
###
###   Pr[X <= x] = exp(-(scale/x)), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvexp <- function (x, scale, log = FALSE)
    .External("do_dpq", "dinvexp", x, scale, log)

##pinvexp <- function(q, scale, tail = TRUE, log = FALSE)
##    .External("do_dpq", "pinvexp", q, scale, tail, log)

##qinvexp <- function(p, scale, tail = TRUE, log = FALSE)
##    .External("do_dpq", "qinvexp", p, scale, tail, log)

rinvexp <- function(n, scale)
    .External("do_random", "rinvexp", n, scale)
