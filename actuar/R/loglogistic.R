### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}loglogistic functions. The loglogistic
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = ((x/scale)^gamma/(1 + (x/scale)^gamma)), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dloglogistic <- function (x, gamma, scale, log)
  .External("do_dpq", "dloglogistic", x, gamma, scale, log)

## ploglogistic <- function(x, gamma, scale)
##     

## qloglogistic <- function(q, gamma, scale)
##    

rloglogistic <- function(n, gamma, scale)
    .External("do_random", "rloglogistic", n, gamma, scale)
