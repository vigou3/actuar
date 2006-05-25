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

diburr <- function (x, shape1, scale, shape2, log = FALSE)
    .External("do_dpq", "diburr", x, shape1, scale, shape2, log = FALSE)    

##piburr <- function(q, shape1, scale, shape2, tail = TRUE, log = FALSE)
##    .External("do_dpq", "piburr", q, shape1, scale, shape2, tail = TRUE, log = FALSE)    

##qiburr <- function(p, shape1, scale, shape2, tail = TRUE, log = FALSE)
##    .External("do_dpq", "qiburr", p, shape1, scale, shape2, tail = TRUE, log = FALSE)

rinverseburr <- function(n, shape1, scale, shape2)
    .External("do_random", "riburr", n, shape1, scale, shape2)
