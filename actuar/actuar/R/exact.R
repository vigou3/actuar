exact <- function(fx, pn)
{
    ## Exact calculation of the total amount of claims probability
    ## function by convolution. Requires a discrete distribution for
    ## the amount of claims.
    ##
    ## ARGUMENTS
    ##
    ## fx: a vector of probabilities of the (discritized) claim amount
    ##     distribution.
    ## pn: a vector of the number of claims probabilities; first element
    ##     must be Pr[N = 0].
    ##
    ## RETURNS
    ##
    ## A vector of probabilities.

    ## Some useful lengths.
    m <- length(fx)      # 1 + maximum claim amount
    n <- length(pn) - 1  # maximum number of claims
    r <- n * m - n + 1   # maximum total amount of claims
    call <- match.call()
    ## Initialization of the output vector.
    fs <- rep(0, r)
    fs[1] <- pn[1]       # Pr[S = 0] = Pr[N = 0]
    ## Convolutions.
    fxc <- 1
    for (i in 1:n)
    {
        pos <- 1:(i * m - i + 1)
        fxc <- convolve(fx, rev(fxc), type="open")
        fs[pos] <- fs[pos] + fxc * pn[i + 1] 
    }
    res <- list(fs = fs, Fs = cumsum(fs), X = 1:length(fs), call = call, FUN = approxfun(cumsum(fs)))
    class(res) <- "AggregateDist"
    res
}
