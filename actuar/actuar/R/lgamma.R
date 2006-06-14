### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}pareto functions. The Loggamma
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = pgamma(log(x), shapelog, 1/ratelog).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dlgamma <- function (x, shapelog, ratelog, log = FALSE)
    .External("do_dpq", "dlgamma", x, shapelog, ratelog, log)

plgamma <- function (q, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "plgamma", q, shapelog, ratelog, lower.tail, log.p)

qlgamma <- function (p, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qlgamma", p, shapelog, ratelog, lower.tail, log.p)

rlgamma <- function(n, shapelog, ratelog)
    .External("do_random", "rlgamma", n, shapelog, ratelog)


