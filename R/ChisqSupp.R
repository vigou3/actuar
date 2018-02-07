### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {m,lev,mgf}chisq functions to compute raw and
### limited moments, and the moment generating function for
### the Chi-square distribution (as defined in R)
###
### See Chapter 17 of Johnson & Kotz, Continuous univariate
### distributions, volume 1, Wiley, 1970
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mchisq <- function(order, df, ncp = 0)
    .External(C_actuar_do_dpq, "mchisq", order, df, ncp, FALSE)

levchisq <- function(limit, df, ncp = 0, order = 1)
    .External(C_actuar_do_dpq, "levchisq", limit, df, ncp, order, FALSE)

mgfchisq <- function(t, df, ncp = 0, log = FALSE)
    .External(C_actuar_do_dpq, "mgfchisq", t, df, ncp, log)
