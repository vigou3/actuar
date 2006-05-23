### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}burr functions. The Burr
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^gamma))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

## dburr <- function (x, shape, scale, gamma)
##     

## pburr <- function(x, shape, scale, gamma)
##     

## qburr <- function(q, shape, scale, gamma)
##     

rburr <- function(n, shape, scale, gamma)
    .External("do_random", "rburr", n, shape, scale, gamma)
