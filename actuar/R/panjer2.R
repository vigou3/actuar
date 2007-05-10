### ===== actuar: an R package for Actuarial Science =====
###
### Panjer recursion formula to compute the approximate aggregate
### claim amount distribution of a portfolio over a period.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair, Louis-Philippe Pouliot and Tommy Ouellet

panjer2 <- function(fx, dist, p0 = NULL, x.scale = 1, ...,
                   TOL = 1e-8, echo = FALSE)
{
    ## Express TOL as a value close to 1.
    TOL <- 1 - TOL

    ## Check whether p0 is a valid probability or not.
    if ( !is.null(p0) ) if ( (p0 < 0) | (p0 > 1) )
        stop("'p0' must be a valid probability (between 0 and 1)")

    ## f_X(0) is no longer needed after the calculation of f_S(0).
    fx0 <- fx[1]
    fx <- fx[-1]

    ## Argument '...' should contain the values of the parameters of
    ## 'dist'.
    par <- list(...)

    ## Distributions are expressed as a member of the (a, b, 0) or (a,
    ## b, 1) families of distributions. Assign parameters 'a' and 'b'
    ## depending of the chosen distribution and compute f_S(0) in
    ## every case, and p1 if p0 is specified in argument.
    if (dist == "geometric")
    {
        dist <- "negative binomial"
        par$size <- 1
    }

    if (dist == "poisson")
    {
        lambda <- par$lambda
        a <- 0
        b <- lambda
        if (is.null(p0))
            fs0 <- exp(-lambda * (1 - fx0))
        else
        {
            fs0 <- p0 + (1 - p0) * (exp(lambda * fx0) - 1)/(exp(lambda) - 1)
            p1 <- (1 - p0) * lambda/(exp(lambda) - 1)
        }
    }
    else if (dist == "negative binomial")
    {
        beta <- 1/(par$prob) - 1
        r <- par$size
        a <- beta/(1 + beta)
        b <- (r - 1) * a
        if (is.null(p0))
            fs0 <- (1 - beta * (fx0 - 1))^(-r)
        else
        {
            fs0 <- p0 + (1 - p0) * ((1 + beta * (1 - fx0))^(-r) - (1 + beta)^(-r))/(1 - (1 + beta)^(-r))
            p1 <- (1 - p0) * r * beta/((1 + beta)^(r + 1) - (1 + beta))
        }
    }
    else if (dist == "binomial")
    {
        m <- par$size
        q <- par$prob
        a <- - q/(1 - q)
        b <- -(m + 1)*a
        if (is.null(p0))
            fs0 <- (1 + q * (fx0 - 1))^m
        else
        {
            fs0 <- p0 + (1 - p0)*((1 + q * (fx0 - 1))^m - (1 - q)^m)/(1 - (1 - q)^m)
            p1 <- (1 - p0) * m * (1 - q)^(m - 1) * q/(1 - (1 - q)^m)
        }
    }
    else if (dist == "logarithmic")
    {
        if (is.null(p0))
            stop("'p0' must be specified with the logarithmic distribution")
        beta <- (1/par$prob) - 1
        a <- beta/(1 + beta)
        b <- -a
        fs0 <- p0 + (1 - p0)*(1 - log(1 - beta(fx0 - 1))/log(1 + beta))
        p1 <- beta/((1 + beta) * log(1 + beta))
    }
    else
        stop("frequency distribution not in the (a, b, 0) or (a, b, 1) families")

    ## If fs0 is equal to zero, the recursion will not start. There is
    ## no provision to automatically cope with this situation in the
    ## current version of this function. Just issue an error message
    ## and let the user do the work by hand.
    if (identical(fs0, 0))
        stop("Pr[S = 0] is numerically equal to 0; impossible to start the recursion")
    if (length(fs0) == 0)
        stop("invalid parameters")

    ## We then use the function .External to do the recursive part of
    ## the Panjer method in C.
    if (is.null(p0)) p1 = 0
    if (is.null(p0)) p0 = -1

    fs <- .External("panjer", p0, p1, fs0, fx0, fx, a, b, TOL, echo)
    fs <- fs[fs != 0]
    
    FUN <- stepfun((0:(length(fs) - 1)) * x.scale, c(0, cumsum(fs)))
    class(FUN) <- c("ecdf", class(FUN))
    assign("fs", fs, env = environment(FUN))
    assign("x.scale", x.scale, env = environment(FUN))
    FUN
}
