### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}transformed beta functions. 
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

## dbetatrans <- function (x, alpha, scale, gamma, tau)
##     

## pbetatrans <- function (x, alpha, scale, gamma, tau)
##     

## qbetatrans <- function (q, alpha, scale, gamma, tau)
##     

rbetatrans <- function (n, alpha, scale, gamma, tau)
    .External("do_random", "rbetatrans", n, alpha, scale, gamma, tau)
