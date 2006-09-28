### ===== actuar: an R package for Actuarial Science =====
###
### Normal and Normal Power Approximation of the total amount of
### claims distribution
###
### See Dayken, Pentikänen and Pesonen, Practical Risk Theory for
### Actuaries, Chapman & Hall, 1994.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

normal <- function(mean, var)
{
    ## Approximate the total amount of claims distribution using the first
    ## two moments.

    call <- match.call()
    FUN <- function(x) pnorm((x - mean)/sqrt(var))
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("call", call, environment(FUN))
    attr(FUN, "source") <- "function(x) pnorm((x - mean)/sqrt(var))"
    comment(FUN) <- "Normal approximation"
    FUN
}

np2 <- function(mean, var, skewness)
{
    ## Approximate the total amount of claims distribution using the first
    ## three moments.
    call <- match.call()

    FUN <- function(x)
    {
        ifelse(x <= mean, NA,
               pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(var) * skewness)) - 3/skewness))                 }

    class(FUN) <- c("aggregateDist", class(FUN))
    comment(FUN) <- "Normal Power approximation"
    assign("call", call, environment(FUN))
    attr(FUN, "source") <- "function(x) pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(var) * skewness)) - 3/skewness))"
    FUN
}
