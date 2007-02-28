### ===== actuar: an R package for Actuarial Science =====
###
### Exact calculation of the aggregate claim amount distribution
### function by convolution. Requires a discrete distribution for
### claim amounts.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Louis-Philippe Pouliot

exact <- function(fx, pn, x.scale = 1)
{
    ## Some useful lengths
    m <- length(fx)                   # 1 + maximum claim amount
    n <- length(pn) - 1               # maximum number of claims
    r <- n * m - n + 1                # maximum total amount of claims

    ## Initialization of the output vector
    fs <- rep(0, r)
    fs[1] <- pn[1]                    # Pr[N = 0]

    ## Convolutions
    fxc <- 1
    for (i in 1:n)
    {
        pos <- seq_len(i * m - i + 1)
        fxc <- convolve(fx, rev(fxc), type = "open")
        fs[pos] <- fs[pos] + fxc * pn[i + 1]
    }

    FUN <- stepfun((0:(length(fs) - 1)) * x.scale, c(0, cumsum(fs)))
    class(FUN) <- c("ecdf", class(FUN))
    assign("fs", fs, environment(FUN))
    assign("x.scale", x.scale, environment(FUN))
    FUN
}
