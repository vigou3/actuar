### ===== actuar: An R Package for Actuarial Science =====
###
###
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>

expint <- function(x)
    .External("actuar_do_dpq", "expint_E1", x, 0, FALSE)


