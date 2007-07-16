### ===== actuar: an R package for Actuarial Science =====
###
### Conditional Tail Expectation for objects of class 'aggregateDist'.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

CTE <- function(...)
    UseMethod("CTE")

CTE.aggregateDist <- function(x, conf.level = 0.99, ...)
{
    label <- comment(x)

    ## The Normal approximation
    if (label == "Normal approximation")
    {
        mean <- get("mean", environment(x))
        variance <- get("variance", environment(x))
        CTEnorm <- exp(-(qnorm(conf.level))^2 / 2) / ( (1 - conf.level) * sqrt(2 * pi) )
        CTE <- CTEnorm * sqrt(variance) + mean
    }

    ## The Normal Power approximation
    else if (label == "Normal Power approximation")
    {
        mean <- get("mean", environment(x))
        variance <- get("variance", environment(x))
        skewness <- get("skewness", environment(x))
        CTEnorm <- exp(-(qnorm(conf.level))^2 / 2) / ( (1 - conf.level) * sqrt(2 * pi) )
        CTE <- ((CTEnorm + 3/skewness)^2 - 9/(skewness^2) - 1) * sqrt(variance) * skewness/6 + mean
    }

    ## The simulation method
    else if (label == "Approximation by simulation")
        CTE = NA
    
    ## The recursive or convolution methods
    else
    {
        pos <- get("x", env = environment(x)) > VaR(x, conf.level)
        CTE <- sum(get("x", env = environment(x))[pos] * get("fs", env = environment(x))[pos]) / (1 - conf.level)
    }
    
    CTE
}
