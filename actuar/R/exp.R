### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}exp functions. The Exponential
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - exp(-x / scale), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mexp <- function(k, scale, log = FALSE)
    .External("do_dpq", "mexp", k, scale, log)

levexp <- function(d, scale, order = 1, log = FALSE)
    .External("do_dpq", "levexp", d, scale, order, log)
