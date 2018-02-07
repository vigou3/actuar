### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {m,lev,mgf}unif functions to compute raw and
### limited moments, and the moment generating function for the
### Uniform distribution (as defined in R).
###
### <http://en.wikipedia.org/wiki/Uniform_distribution_%28continuous%29>
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

munif <- function(order, min = 0, max = 1)
    .External(C_actuar_do_dpq, "munif", order, min, max, FALSE)

levunif <- function(limit, min = 0, max =1, order = 1)
    .External(C_actuar_do_dpq, "levunif", limit, min, max, order, FALSE)

mgfunif <- function(t, min = 0, max = 1, log = FALSE)
    .External(C_actuar_do_dpq, "mgfunif", t, min, max, log)
