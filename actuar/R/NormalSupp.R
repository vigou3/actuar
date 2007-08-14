### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}exp functions to compute raw and limited
### moments, and the moment generating function
### for the Normal distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter 13 of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mnorm <- function(order, mean = 0, sd = 1)
    .External("do_dpq", "mnorm", order, mean, sd, FALSE)

mgfnorm <- function(t, mean = 0, sd = 1, log = FALSE)
    .External("do_dpq", "mgfnorm", t, mean, sd, log)
