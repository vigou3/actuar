### ===== actuar: an R package for Actuarial Science =====
###
### Value at Risk for objects of class 'aggregateDist'.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

VaR <- function(x, ...)
    UseMethod("VaR")

VaR.aggregateDist <- function(x, conf.level = 0.99, ...)
    quantile(x, conf.level)
