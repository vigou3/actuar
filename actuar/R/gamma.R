### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}gamma functions.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mgamma <- function(order, shape, rate = 1, scale = 1/rate)
    .External("do_dpq", "mgamma", order, shape, scale, FALSE)

levgamma <- function(limit, shape, rate = 1, scale = 1/rate, order = 1)
    .External("do_dpq", "levgamma", limit, shape, scale, order, FALSE)
