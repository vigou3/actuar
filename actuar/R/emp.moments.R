### ===== actuar: an R package for Actuarial Science =====
###
### Empirical moments for individual and grouped data.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

emp.moments <- function(x, order = 1, ...)
  UseMethod("emp.moments")

emp.moments.default <- function(x, order = 1, ...)
{
    if (any(order < 0))
      stop("order must be positive")
    
    sapply(order, function(n) mean(x^n))
}

emp.moments.grouped.data <- function(x, order = 1, ...)
{
    if (any(!is.finite(x$cj)))
      stop("infinite value")
    if (any(order < 0))
      stop("order must be positive")
    
    mom.1 <- function(k, x)
    {
        n <- sum(x$nj)
        sum(diff(x$cj^(k+1)) / diff(x$cj) / (k+1) * x$nj[-1] / n)
    }
    sapply(order, mom.1, x = x)
}

