### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse exponential functions. The inverse exponential
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = exp(-(scale/x)), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinverseexp <- function (x, scale, log)
  .External("do_dpq", "dinverseexp", x, scale, log)    

## pinverseexp <- function(x, scale)
##    

## qinverseexp <- function(q, scale)
##  

rinverseexp <- function(n, scale)
    .External("do_random", "rinverseexp", n, scale)
