### ===== actuar: an R package for Actuarial Science =====
###
### Exact calculation of the total amount of claims distribution
### function by convolution. Requires a discrete distribution for
### the amount of claims.
###
### ARGUMENTS
###
### fx: a vector of probabilities of the (discritized) claim amount
###     distribution.
### pn: a vector of the number of claims probabilities; first element
###     must be Pr[N = 0].
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### and Louis-Philippe Pouliot

exact <- function(x.scale = 1, fx, pn)
{
    ## Some useful lengths.
    call <- match.call()
    m <- length(fx)      # 1 + maximum claim amount
    n <- length(pn) - 1  # maximum number of claims
    r <- n * m - n + 1   # maximum total amount of claims
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
    FUN <- stepfun((0:(length(fs) - 1)) * x.scale, c(0, cumsum(fs)))
    class(FUN) <- c("aggregateDist", "ecdf", class(FUN))
    assign("call", call, environment(FUN))
    assign("fs", fs, environment(FUN))
    assign("x.scale", x.scale, environment(FUN))
    comment(FUN) <- "Direct calculation"
    FUN
}

