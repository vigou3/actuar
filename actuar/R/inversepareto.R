
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

dipareto <- function (x, shape, scale, log = FALSE)
     .External("do_dpq", "dipareto", x, shape, scale, log)     

##pipareto <- function(q, shape, scale, tail = TRUE, log = FALSE)
##    .External("do_dpq", "pipareto", q, shape, scale, tail, log) 

##qipareto <- function(p, shape, scale, tail = TRUE, log = FALSE)
##     .External("do_dpq", "qipareto", p, shape, scale, tail, log) 

ripareto <- function(n, shape, scale)
     .External("do_random", "ripareto", n, shape, scale)
