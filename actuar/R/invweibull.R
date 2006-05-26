### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse weibull functions. The inverse Weibull
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - exp(-(x/scale)^tau), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvweibull <- function (x, scale, shape, log = FALSE)
    .External("do_dpq", "dinvweibull", x, scale, shape, log)

#pinvweibull <- function(q, scale, shape, tail = TRUE, log = FALSE)
#     .External("do_dpq", "pinvweibull", q, scale, shape, tail, log)

#qinvweibull <- function(p, scale, shape, tail = TRUE, log = FALSE)
#     .External("do_dpq", "qinvweibull", p, scale, shape, tail, log)

rinvweibull <- function(n, scale, shape)
    .External("do_random", "rinvweibull", n, scale, shape)

## Aliases
dlgompertz <- dinvweibull
#plgompertz <- pinvweibull
#qlgompertz <- qinvweibull
rlgompertz <- rinvweibull
