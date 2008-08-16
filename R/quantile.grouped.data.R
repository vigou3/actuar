### ===== actuar: an R package for Actuarial Science =====
###
### Quantiles (inverse of the ogive) for grouped data
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

quantile.grouped.data <- function(x, probs = seq(0, 1, 0.25),
                                  names = TRUE, ...)
{
    y <- x[, 2]
    x <- eval(expression(cj), env = environment(x))

    ## Inverse of the ogive
    fun <- approxfun(c(0, cumsum(y))/sum(y), x,
                     yleft = min(x), yright = max(x),
                     method = "linear", ties = "ordered")

    ## Quantiles
    res <- fun(probs)

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", wid = 1, digits = dig)
    }
    res
}
