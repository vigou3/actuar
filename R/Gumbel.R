### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}gumbel functions to compute
### characteristics of the Gumbel distribution. The version used in
### these functions has cumulative distribution function
###
###   Pr[X <= x] = exp(-exp(-(x - alpha)/beta)), -Inf < x < Inf.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dgumbel <- function (x, alpha, beta, log = FALSE)
    .External("actuar_do_dpq", "dgumbel", x, alpha, beta, log)

pgumbel <- function(q, alpha, beta, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pgumbel", q, alpha, beta, lower.tail, log.p)

qgumbel <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qgumbel", p, alpha, beta, lower.tail, log.p)

rgumbel <- function(n, alpha, beta)
    .External("actuar_do_random", "rgumbel", n, alpha, beta)

mgumbel <- function(order, alpha, beta)
    .External("actuar_do_dpq", "mgumbel", order, alpha, beta, FALSE)

mgfgumbel <- function(t, alpha, beta, log = FALSE)
    .External("actuar_do_dpq", "mgfgumbel", t, alpha, beta, log)
