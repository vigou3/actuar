### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}weibull functions. 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mweibull <- function(k, scale, shape, log = FALSE)
    .External("do_dpq", "mweibull", k, scale, shape, log)

levweibull <- function(d, scale, shape, order = 1, log = FALSE)
    .External("do_dpq", "levweibull", d, scale, shape, order, log)
