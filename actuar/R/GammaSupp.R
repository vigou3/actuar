### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}gamma functions to compute raw and
### limited moments, and the moment generating function for 
### the Gamma distribution (as defined in R) 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter 17 of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
### AUTHORS:  Mathieu Pigeon, Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

mgamma <- function(order, shape, rate = 1, scale = 1/rate)
    .External("do_dpq", "mgamma", order, shape, scale, FALSE)

levgamma <- function(limit, shape, rate = 1, scale = 1/rate, order = 1)
    .External("do_dpq", "levgamma", limit, shape, scale, order, FALSE)

mgfgamma <- function(x, shape, rate = 1, scale = 1/rate, log= FALSE)
    .External("do_dpq", "mgfgamma", x, shape, scale, log)
