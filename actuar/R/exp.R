### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}exp functions to compute raw and limited
### moments for the Exponential distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mexp <- function(order, rate = 1)
    .External("do_dpq", "mexp", order, 1/rate, FALSE)

levexp <- function(limit, rate = 1, order = 1)
    .External("do_dpq", "levexp", limit, 1/rate, order, FALSE)
