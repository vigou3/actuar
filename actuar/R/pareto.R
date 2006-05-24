### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}pareto functions. The Pareto
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (scale/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto <- function (x, shape, scale, log)
    .External("do_dpq", "dpareto", x, shape, scale, log)
     
ppareto <- function (x, shape, scale, tail, log)
    .External("do_dpq", "ppareto", x, shape, scale, tail, log)

qpareto <- function (q, shape, scale, tail, log)
    .External("do_dpq", "qpareto", q, shape, scale, tail, log)

rpareto <- function(n, shape, scale)
    .External("do_random", "rpareto", n, shape, scale)
