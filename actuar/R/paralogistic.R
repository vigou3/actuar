### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}paralogistic functions. The paralogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^alpha))^alpha, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dparalogistic <- function (x, alpha, scale, log)
  .External("do_dpq", "dparalogistic", x, alpha, scale, log)

## pparalogistic <- function(x, alpha, scale)
##    

## qparalogistic <- function(q, alpha, scale)
##

rparalogistic <- function(n, alpha, scale)
    .External("do_random", "rparalogistic", n, alpha, scale)
