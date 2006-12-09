### ===== actuar: an R package for Actuarial Science =====
###
### 'mean' method for 'aggregateDist' objects
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot
mean.aggregateDist <- function(x, ...)
{
    ## Simply return the value of the true mean given in argument in
    ## the case of the Normal and Normal Power approximations.
    if (comment(x) %in%
        c("Normal approximation", "Normal Power approximation"))
        return(eval(expression(mean), environment(x)))

    ## For the recursive, exact and simulation methods, compute the
    ## mean from the stepwise cdf. For the first two, use the pmf
    ## saved in the environment of the object.
    val <- get("x", environment(x))
    prob <-
        if (comment(x) == "Approximation by simulation")
            c(x[1], diff(get("y", environment(x))))
        else
            get("fs", environment(x))
    drop(crossprod(val, prob))
}
