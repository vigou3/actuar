### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}gamma functions to compute raw and
### limited moments, and the moment generating function for 
### the Gamma distribution (as defined in package SuppDists) 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
### See Wikipedia, Inverse Gaussian for the "mgf"
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

minvGauss <- function(order, nu, lambda )
    .External("do_dpq", "minvGauss", order, nu, lambda, FALSE)

levinvGauss <- function(limit, nu, lambda, order = 1)
    .External("do_dpq", "levinvGauss", limit,  nu, lambda, order, FALSE)

mgfinvGauss <- function(t, nu, lambda, log= FALSE)
    .External("do_dpq", "mgfinvGauss", t, nu, lambda, log)
