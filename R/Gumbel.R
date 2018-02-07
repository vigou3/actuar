### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r,m,lev}gumbel functions to compute
### characteristics of the Gumbel distribution. The version used in
### these functions has cumulative distribution function
###
###   Pr[X <= x] = exp(-exp(-(x - alpha)/scale)), -Inf < x < Inf.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dgumbel <- function (x, alpha, scale, log = FALSE)
    .External(C_actuar_do_dpq, "dgumbel", x, alpha, scale, log)

pgumbel <- function(q, alpha, scale, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pgumbel", q, alpha, scale, lower.tail, log.p)

qgumbel <- function(p, alpha, scale, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qgumbel", p, alpha, scale, lower.tail, log.p)

rgumbel <- function(n, alpha, scale)
    .External(C_actuar_do_random, "rgumbel", n, alpha, scale)

mgumbel <- function(order, alpha, scale)
    .External(C_actuar_do_dpq, "mgumbel", order, alpha, scale, FALSE)

mgfgumbel <- function(t, alpha, scale, log = FALSE)
    .External(C_actuar_do_dpq, "mgfgumbel", t, alpha, scale, log)
