### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse pareto functions. The inverse Pareto
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = (x/(x + scale))^tau, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

## dinversepareto <- function (x, tau, scale)
##     

## pinversepareto <- function(x, tau, scale)
##     

## qinversepareto <- function(q, tau, scale)
##    

rinversepareto <- function(n, tau, scale)
    .External("do_random", "rinversepareto", n, tau, scale)
