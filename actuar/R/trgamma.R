### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}transformed gamma functions. The transformed gamma
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = pgamma((x / scale)^shape2 , shape1, scale), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

#dtrgamma <- function (x, shape1, scale, shape2, log = FALSE)
#    .External("do_dpq", "dtrgamma", x, shape1, scale, shape2, log)

#ptrgamma <- function(q, shape1, scale, shape2, tail = TRUE, log = FALSE)
#    .External("do_dpq", "ptrgamma", q, shape1, scale, shape2, tail, log)     

#qtrgamma <- function(p, shape1, scale, shape2, tail = TRUE, log = FALSE)
#    .External("do_dpq", "qtrgamma", p, shape1, scale, shape2, tail, log)   

#rtrgamma <- function(n, shape1, scale, shape2)
#    .External("do_random", "rtrgamma", n, shape1, scale, shape2)
