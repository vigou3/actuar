### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}gamma functions. 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mgamma <- function(k, shape, scale, log = FALSE)
    .External("do_dpq", "mgamma", k, shape, scale, log)

levgamma <- function(d, shape, scale, order = 1, log = FALSE)
    .External("do_dpq", "levgamma", d, shape, scale, order, log)
