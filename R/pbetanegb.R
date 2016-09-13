### ===== actuar: An R Package for Actuarial Science =====
###
### Function to compute the incomplete beta function when b < 0, b !=
### -1, -2, ... and a > 1+ floor(-b). Used mostly at the C level to
### compute the limited expected value for distributions of the
### transformed beta family.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Wiley.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

pbetanegb <- function (x, a, b)
    .External("actuar_do_dpq", "pbetanegb", x, a, b, FALSE)


