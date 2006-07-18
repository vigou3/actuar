### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}inverse weibull functions. The inverse Weibull
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - exp(-(x/scale)^shape), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvweibull <- function (x, scale, shape, log = FALSE)
    .External("do_dpq", "dinvweibull", x, scale, shape, log)

pinvweibull <- function(q, scale, shape, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "pinvweibull", q, scale, shape, lower.tail, log.p)

qinvweibull <- function(p, scale, shape, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qinvweibull", p, scale, shape, lower.tail, log.p)

rinvweibull <- function(n, scale, shape)
    .External("do_random", "rinvweibull", n, scale, shape)

minvweibull <- function(k, scale, shape, log = FALSE)
     .External("do_dpq", "minvgamma", k, scale, shape, log)

levinvweibull <- function(d, scale, shape, order = 1, log = FALSE)
     .External("do_dpq", "levinvweibull", d, scale, shape, order, log)

## Aliases
dlgompertz <- dinvweibull
plgompertz <- pinvweibull
qlgompertz <- qinvweibull
rlgompertz <- rinvweibull
mlgompertz <- minvweibull
levlgompertz <- levinvweibull
