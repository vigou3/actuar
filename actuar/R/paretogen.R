### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}generalized pareto functions. 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dparetogen <- function (x, shape, scale, tau, log)
  .External("do_dpq", "dparetogen", x, shape, scale, tau, log)

## pparetogen <- function(x, shape, scale, tau)
##     

## qparetogen <- function(q, shape, scale, tau)
##     

rparetogen <- function(n, shape, scale, tau)
    .External("do_random", "rparetogen", n, shape, scale, tau)
