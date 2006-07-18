### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}pareto functions. The Pareto
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

ppareto <- function (q, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "ppareto", q, shape, scale, lower.tail, log.p)

qpareto <- function (p, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qpareto", p, shape, scale, lower.tail, log.p)

rpareto <- function(n, shape, scale)
    .External("do_random", "rpareto", n, shape, scale)

mpareto <- function(k, shape, scale, log = FALSE)
     .External("do_dpq", "mpareto", k, shape, scale, log)

levpareto <- function(d, shape, scale, order = 1, log = FALSE)
     .External("do_dpq", "levpareto", d, shape, scale, order, log)

## Aliases
dpareto2 <- dpareto
ppareto2 <- ppareto
qpareto2 <- qpareto
rpareto2 <- rpareto
mpareto2 <- mpareto
levpareto2 <- levpareto
