### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}burr functions. The Burr
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^shape2))^shape1, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dburr <- function (x, shape1, scale, shape2, log = FALSE)
    .External("do_dpq", "dburr", x, shape1, scale, shape2, log)

pburr <- function(q, shape1, scale, shape2, tail = TRUE, log = FALSE)
    .External("do_dpq", "pburr", q, shape1, scale, shape2, tail, log)     

qburr <- function(p, shape1, scale, shape2, tail = TRUE, log = FALSE)
    .External("do_dpq", "qburr", p, shape1, scale, shape2, tail, log)   

rburr <- function(n, shape1, scale, shape2)
    .External("do_random", "rburr", n, shape1, scale, shape2)
