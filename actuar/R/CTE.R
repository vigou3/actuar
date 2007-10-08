### ===== actuar: an R package for Actuarial Science =====
###
### Conditional Tail Expectation for objects of class 'aggregateDist'.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

CTE <- function(x, ...)
    UseMethod("CTE")

CTE.aggregateDist <- function(x, conf.level = c(0.9, 0.95, 0.99),
                              names = TRUE, ...)
{
    label <- comment(x)

    ## Normal approximation; an exact formula is available
    if (label == "Normal approximation")
    {
        m <- get("mean", environment(x))
        sd <- sqrt(get("variance", environment(x)))
        res <- m + sd * exp(-(qnorm(conf.level))^2 / 2) /
            ((1 - conf.level) * sqrt(2 * pi))
    }
    ## Normal Power approximation; no explicit formula so revert to
    ## numerical integration.
    else if (label == "Normal Power approximation")
    {
        m <- get("mean", envir = environment(x))
        sd <- sqrt(get("variance", envir = environment(x)))
        sk <- get("skewness", envir = environment(x))

        f <- function(x)
        {
            y <- sqrt(1 + 9/sk^2 + 6 * (x - m)/(sd * sk))
            3 * x * dnorm(y - 3/sk) / (sd * sk * y)
        }

        res <- sapply(quantile(x, conf.level),
                      function(x) integrate(f, x, Inf)$value) /
                          (1 - conf.level)
    }
    ## Recursive method, simulation and convolutions; each yield a
    ## step function that can be used to make calculations.
    else
    {
        val <- get("x", envir = environment(x))
        prob <- get("fs", envir = environment(x))
        f <- function(a)
        {
            pos <- val > VaR(x, a)
            drop(crossprod(val[pos], prob[pos])) / (1 - a)
        }
        res <- sapply(conf.level, f)
    }

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * conf.level, "%", sep = ""),
                              format = "fg", wid = 1, digits = dig)
    }
    res
}
