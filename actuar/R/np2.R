np2 <- function(mean, var, skewness)
{

    ## Evaluate a distribution function using the Normal Power
    ## approximation.
    ##
    ## ARGUMENTS
    ##
    ##        x: value at which the function is evaluated
    ##     mean: true mean of the distribution
    ##      var: true variance of the distribution
    ## skewness: true skewness of the distribution
    ##
    ## RETURNS
    ##
    ## 'NA' if x <= mean, probability otherwise.
    call <- match.call()
    
    FUN <- function(x)
    {
        ifelse(x <= mean, NA,
               pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(var) * skewness)) - 3/skewness))                 }
    
    class(FUN) <- c("aggregateDist", class(FUN))
    #assign("label", "Normal Power approximation", environment(FUN))
    comment(FUN) <- "Normal Power approximation"
    assign("call", call, environment(FUN))
    attr(FUN, "source") <- "function(x) pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(var) * skewness)) - 3/skewness))"
    FUN    
}




