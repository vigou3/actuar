### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}single-parameter pareto functions. The single-parameter ### Pareto distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (scale/x)^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dspareto <- function (x, shape, scale, log)
    .External("do_dpq", "dspareto", x, shape, scale, log)

## pspareto <- function(x, shape, scale)
##    

## qspareto <- function(q, shape, scale)
##   

rspareto <- function(n, shape, scale)
    .External("do_random", "rspareto", n, shape, scale)
