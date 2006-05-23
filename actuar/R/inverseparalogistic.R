### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse paralogistic functions. The inverse paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^tau/(1 + (x/scale)^tau))^tau, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

## dinverseparalogistic <- function (x, tau, scale)
##    

## pinverseparalogistic <- function(x, tau, scale)
##     

## qinverseparalogistic <- function(q, tau, scale)
##     

rinverseparalogistic <- function(n, tau, scale)
    .External("do_random", "rinverseparalogistic", n, tau, scale)
