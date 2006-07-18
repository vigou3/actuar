### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}lnorm functions.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mlnorm <- function(k, mu, sigma, log = FALSE)
    .External("do_dpq", "mlnorm", k, mu, sigma, log)

levlnorm <- function(d, mu, sigma, order = 1, log = FALSE)
    .External("do_dpq", "levlnorm", d, mu, sigma, order, log)
