### ===== actuar: an R package for Actuarial Science =====
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
{
    if (all(as.integer(order) == order))
        .External("do_dpq", "mchisq", order, df, ncp, FALSE)
    else
        stop("'order' must be integers")
}

levchisq <- function(limit, df, ncp = 0, order = 1)
{
    if (all(ncp == 0))
        .External("do_dpq", "levchisq", limit, df, ncp, order, FALSE)
    else
        stop("non central chi-square not supported")
}

mgfchisq <- function(x, df, ncp = 0, log = FALSE)
    .External("do_dpq", "mgfchisq", x, df, ncp, log)
