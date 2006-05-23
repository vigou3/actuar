### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse burr functions. The inverse Burr
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^gamma/(1 + (x/scale)^gamma))^tau x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

## dinverseburr <- function (x, tau, scale, gamma)
##    

## pinverseburr <- function(x, tau, scale, gamma)
##     

## qinverseburr <- function(q, tau, scale, gamma)
##     

rinverseburr <- function(n, tau, scale, gamma)
    .External("do_random", "rinverseburr", n, tau, scale, gamma)
