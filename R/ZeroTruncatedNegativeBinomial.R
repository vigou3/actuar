### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}ztnbinom functions to compute
### characteristics of the Zero Truncated Negative Binomial
### distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dztnbinom <- function (x, size, prob, log = FALSE)
    .External(C_actuar_do_dpq, "dztnbinom", x, size, prob, log)

pztnbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pztnbinom", q, size, prob, lower.tail, log.p)

qztnbinom <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qztnbinom", p, size, prob, lower.tail, log.p)

rztnbinom <- function(n, size, prob)
    .External(C_actuar_do_random, "rztnbinom", n, size, prob)

## not exported; for internal use in panjer()
pgfztnbinom <- function(x, size, prob)
    expm1(-size * log1p(x * (prob - 1)))/expm1(-size * log(prob))
