### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}invgauss functions to compute raw and
### limited moments, and the moment generating function for
### the Inverse Gaussian distribution (as defined in package SuppDists)
###
### See Part 2 of Chhikara and Folks, The inverse Gaussian
### distribution: theory, methodology and application, Decker, 1989;
### Part 1 of Seshadri, The inverse Gaussian distribution: statistical
### theory and applications, Springer, 1989
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

minvGauss <- function(order, nu, lambda)
{
    if (all(as.integer(order) == order))
        .External("do_dpq", "minvGauss", order, nu, lambda, FALSE)
    else
        stop("'order' must be integers")
}

levinvGauss <- function(limit, nu, lambda, order = 1)
{
    if (order == 1)
        .External("do_dpq", "levinvGauss", limit,  nu, lambda, order, FALSE)
    else
        stop("orders other than 1 are not supported")
}

mgfinvGauss <- function(x, nu, lambda, log= FALSE)
    .External("do_dpq", "mgfinvGauss", x, nu, lambda, log)

## Aliases
minvgauss <- minvGauss
levinvgauss <- levinvGauss
mgfinvgauss <- mgfinvGauss
