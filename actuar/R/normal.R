normal <- function(mean, var)
{
    ## Evaluate a distribution function using the normal
    ## approximation.
    ## Takes two arguments: the true mean and the true variance
    
    call <- match.call()
    FUN <- function(x) pnorm((x - mean)/sqrt(var))
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("call", call, environment(FUN))
    #assign("label", "Normal approximation", environment(FUN))
    attr(FUN, "source") <- "function(x) pnorm((x - mean)/sqrt(var))"
    comment(FUN) <- "Normal approximation"
    FUN
}

