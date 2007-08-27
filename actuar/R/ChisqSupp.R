### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}chisq functions to compute raw and
### limited moments, and the moment generating function for 
### the Chi-square distribution (as defined in R) 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Chapter 17 of Johnson & Kotz, Loss Distributions, Wiley, 1970
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mchisq <- function(order, df, ncp = 0)
{
    if( as.integer(order) == order )
        .External("do_dpq", "mchisq", order, df, ncp, FALSE)
    else
        stop("non integer order is not supported")
}

levchisq <- function(limit, df, ncp = 0, order = 1)
{
    if( ncp == 0 )
        .External("do_dpq", "levchisq", limit, df, ncp, order, FALSE)
    else
        stop("non central (ncp>0) chi-square is not supported")
}
mgfchisq <- function(x, df, ncp = 0, log= FALSE)
    .External("do_dpq", "mgfchisq", x, df, ncp, log)
