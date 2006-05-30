### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}pareto functions. The Pareto
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (scale/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto <- function (x, shape, scale, log = FALSE)
    .External("do_dpq", "dpareto", x, shape, scale, log)

ppareto <- function (q, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "ppareto", q, shape, scale, tail, log)

qpareto <- function (p, shape, scale, tail = TRUE, log = FALSE)
    .External("do_dpq", "qpareto", p, shape, scale, tail, log)

rpareto <- function(n, shape, scale)
    .External("do_random", "rpareto", n, shape, scale)

## Aliases
dpareto2 <- dpareto
ppareto2 <- ppareto
qpareto2 <- qpareto
rpareto2 <- rpareto
