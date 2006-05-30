### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}transformed beta functions.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dtrbeta <- function (x, shape1, scale, shape2, shape3, log = FALSE)
    .External("do_dpq", "dtrbeta", x, shape1, scale, shape2, shape3, log)

ptrbeta <- function (q, shape1, scale, shape2, shape3, tail = TRUE, log = FALSE)
    .External("do_dpq", "ptrbeta", q, shape1, scale, shape2, shape3, tail, log)

qtrbeta <- function (p, shape1, scale, shape2, shape3, tail = TRUE, log = FALSE)
    .External("do_dpq", "qtrbeta", p, shape1, scale, shape2, shape3, tail, log)

rtrbeta <- function (n, shape1, scale, shape2, shape3)
    .External("do_random", "rtrbeta", n, shape1, scale, shape2, shape3)

## Aliases
dpearson4 <- dtrbeta
ppearson4 <- ptrbeta
qpearson4 <- qtrbeta
rpearson4 <- rtrbeta
