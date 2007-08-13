### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}exp functions to compute raw and limited
### moments,  and the moment generating function
### for the Exponential distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter 18 of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
###
### AUTHORS:  Mathieu Pigeon, Christophe Dutang
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

mexp <- function(order, rate = 1)
    .External("do_dpq", "mexp", order, 1/rate, FALSE)

levexp <- function(limit, rate = 1, order = 1)
    .External("do_dpq", "levexp", limit, 1/rate, order, FALSE)

mgfexp <- function(x, rate = 1, log = FALSE)
    .External("do_dpq", "mgfexp", x, 1/rate, log)
