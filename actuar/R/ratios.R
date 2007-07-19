### ===== actuar: an R package for Actuarial Science =====
###
### Evaluate the ratios of the aggregate claim amounts to
### natural weights for objects of class "simpf".
###
### AUTHORS: Tommy Ouellet,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

ratios <- function(object, ...)
    UseMethod("ratios")
