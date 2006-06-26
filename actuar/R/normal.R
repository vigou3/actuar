normal <- function(mean, var)
{
    FUN <- function(x) pnorm((x - eval(mean))/sqrt(eval(var)))
    class(FUN) <- c("aggregateDist", class(FUN))
    FUN
}

