### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {m,lev}beta functions to compute raw and limited
### moments for the Beta distribution (as defined in R). The
### noncentral beta distribution is _not_ supported.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

mbeta <- function(order, shape1, shape2)
    .External(C_actuar_do_dpq, "mbeta", order, shape1, shape2, FALSE)

levbeta <- function(limit, shape1, shape2, order = 1)
    .External(C_actuar_do_dpq, "levbeta", limit, shape1, shape2, order, FALSE)
