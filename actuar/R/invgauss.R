### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r}inverse Gaussian functions.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvgauss <- function (x, mean, scale, log = FALSE)
     .External("do_dpq", "dinvgauss", x, mean, scale, log)

pinvgauss <- function(q, mean, scale, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpq", "pinvgauss", q, mean, scale, lower.tail, log.p)

qinvgauss <- function(p, mean, scale, lower.tail = TRUE, log.p = FALSE)
     .External("do_dpq", "qinvgauss", p, mean, scale, lower.tail, log.p)

rinvgauss <- function(n, mean, scale)
     .External("do_random", "rinvgauss", n, mean, scale)
