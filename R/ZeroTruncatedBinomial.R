### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Definition of the {d,p,q,r}ztbinom functions to compute
### characteristics of the Zero Truncated Binomial distribution.
###
### See Appendix B of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dztbinom <- function (x, size, prob, log = FALSE)
    .External(C_actuar_do_dpq, "dztbinom", x, size, prob, log)

pztbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "pztbinom", q, size, prob, lower.tail, log.p)

qztbinom <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE)
    .External(C_actuar_do_dpq, "qztbinom", p, size, prob, lower.tail, log.p)

rztbinom <- function(n, size, prob)
    .External(C_actuar_do_random, "rztbinom", n, size, prob)

## not exported; for internal use in panjer()
pgfztbinom <- function(x, size, prob)
{
    qn <- (1 - prob)^size
    (exp(size * log1p(prob * (x - 1))) - qn)/(1 - qn)
}
