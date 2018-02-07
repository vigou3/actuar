### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {m,lev}lnorm functions.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mlnorm <- function(order, meanlog = 0, sdlog = 1)
    .External(C_actuar_do_dpq, "mlnorm", order, meanlog, sdlog, FALSE)

levlnorm <- function(limit, meanlog = 0, sdlog = 1, order = 1)
    .External(C_actuar_do_dpq, "levlnorm", limit, meanlog, sdlog, order, FALSE)
