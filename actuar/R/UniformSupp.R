### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}exp functions to compute raw and limited
### moments, and the moment generating function
### for the Uniform distribution (as defined in R).
###
### Wikipedia, The uniform distribution (continuous)
###
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

munif <- function(order, min = 0, max = 1)
    .External("do_dpq", "munif", order, min, max, FALSE)

levunif <- function(limit, min = 0, max =1, order = 1)
    .External("do_dpq", "levunif", limit, min, max, order, FALSE)

mgfunif <- function(t, min = 0, max = 1, log = FALSE)
    .External("do_dpq", "mgfunif", t, min, max, log)
