### Definition of the {d,p,q,r}pareto functions. The Pareto
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (scale/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.

## dpareto <- function (x, alpha, lambda)
##     .C("R_dpareto", as.double(x), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

## ppareto <- function(x, alpha, lambda)
##     .C("R_ppareto", as.double(x), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

## qpareto <- function(q, alpha, lambda)
##     .C("R_qpareto", as.double(q), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

rpareto <- function(n, shape, scale)
    .C("rpareto_sym", as.integer(n), as.double(shape), as.double(scale),
       y = double(n))$y
