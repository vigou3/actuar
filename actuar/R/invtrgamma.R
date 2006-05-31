### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse transformed gamma functions. The Burr
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^gamma))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvtrgamma <- function (x, shape1, scale, shape2, log = FALSE)
    .External("do_dpq", "dinvtrgamma", x, shape1, scale, shape2, log)

pinvtrgamma <- function(q, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "pinvtrgamma", q, shape1, scale, shape2, lower.tail, log.p)     

qinvtrgamma <- function(p, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "qinvtrgamma", p, shape1, scale, shape2, lower.tail, log.p)   

rinvtrgamma <- function(n, shape1, scale, shape2)
    .External("do_random", "rinvtrgamma", n, shape1, scale, shape2)
