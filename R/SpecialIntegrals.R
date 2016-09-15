### ===== actuar: An R Package for Actuarial Science =====
###
### Function to compute the incomplete beta function when b < 0, b !=
### -1, -2, ... and a > 1+ floor(-b). Used mostly at the C level to
### compute the limited expected value for distributions of the
### transformed beta family.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

betaint <- function(x, a, b)
    .External("actuar_do_dpq", "betaint", x, a, b, FALSE)

expint <- function(x)
    .External("actuar_do_dpq", "expint", x, FALSE)

gammaint <- function(x, a)
    .External("actuar_do_dpq", "gammaint", x, a, FALSE)

