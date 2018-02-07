### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,r}zmlogarithmic functions to compute
### characteristics of the zero modified logarithmic distribution. See
### ./Logarithmic.R for details on the parametrization.
###
### See p. 93 of Klugman, Panjer & Willmot, Loss Models, Fourth
### Edition, Wiley, 2012.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dzmlogarithmic <- function(x, prob, p0, log = FALSE)
    .External(C_actuar_do_dpq, "dzmlogarithmic", x, prob, p0, log)

pzmlogarithmic <- function(q, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pzmlogarithmic", q, prob, p0, lower.tail, log.p)

qzmlogarithmic <- function(p, prob, p0, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qzmlogarithmic", p, prob, p0, lower.tail, log.p)

rzmlogarithmic <- function(n, prob, p0)
    .External(C_actuar_do_random, "rzmlogarithmic", n, prob, p0)
