### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse weibull functions. The inverse Weibull
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - exp(-(x/scale)^tau), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dlgompertz <- function (x, scale, shape, log = FALSE)
     .External("do_dpq", "dlgompertz", x, scale, shape, log)     

##plgompertz <- function(q, scale, shape, tail = TRUE, log = FALSE)
##     .External("do_dpq", "plgompertz", q, scale, shape, tail, log) 

##qlgompertz <- function(p, scale, shape, tail = TRUE, log = FALSE)
##     .External("do_dpq", "qlgompertz", p, scale, shape, tail, log) 

rlgompertz <- function(n, scale, shape)
     .External("do_random", "rlgompertz", n, scale, shape)
