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

## dinverseweibull <- function (x, scale, tau)
##     

## pinverseweibull <- function(x, scale, tau)
##     

## qinverseweibull <- function(q, scale, tau)
##    

rinverseweibull <- function(n, scale, tau)
    .External("do_random", "rinverseweibull", n, scale, tau)
