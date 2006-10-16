### ===== actuar: an R package for Actuarial Science =====
###
### Empirical moments for individual and grouped data.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

emoments <- function(x, order = 1, ...)
  UseMethod("emp.moments")

emoments.default <- function(x, order = 1, ...)
{
    if (any(order < 0))
      stop("'order' must be non negative")

    sapply(order, function(n) mean(x^n))
}

emoments.grouped.data <- function(x, order = 1, ...)
{
    ## Sanity check
    if (!inherits(x, "grouped.data"))
        stop("wrong method")

    if (any(order < 0))
        stop("'order' must be non negative")

    cj <- eval(expression(cj), env = environment(x))

    if (any(!is.finite(cj)))
        stop("cannot compute moments with an infinite group")

    cjt <- cj[-length(cj)]

    sapply(x[-1], function(x) drop(crossprod(x, midpoints)))


    fun <- fxunction(x, k)
    {
        n <- sum(x)
        sum(x * diff(cj^(k+1)) / diff(x$cj) / (k+1) * x$nj[-1] / n)
    }
    sapply(order, mom.1, x = x)
}
