np2 <- function(x, mean, var, skewness)
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
    Fs <- ifelse(x <= mean, NA,
                 pnorm(sqrt(1 + 9/skewness^2 + 6 * (x - mean)/(sqrt(var) * skewness)) - 3/skewness))
    res <- list(fs = c(NA, diff(Fs)), Fs = Fs, X = 1:length(Fs), call = call, FUN = approxfun(Fs))
    class(res) <- "AggregateDist"
    res    
}
