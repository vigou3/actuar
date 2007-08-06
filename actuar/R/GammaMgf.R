### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the moment generating function for the 
### Gamma distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>


mgfgamma <- function(t, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("do_dpq", "mgfgamma", t, shape, scale, order, log)
