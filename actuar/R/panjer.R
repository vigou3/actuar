### ===== actuar: an R package for Actuarial Science =====
###
### Panjer recursion formula to compute the approximate aggreagte
### claim amount distribution of a portfolio over a period.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair and Louis-Philippe Pouliot

panjer <- function(fx, model.freq, p0 = NULL, x.scale = 1, TOL = 1e-8, echo = FALSE)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()

    ## Express TOL as a value close to 1.
    TOL <- 1 - TOL

    ## f_X(0) is no longer needed after the calculation of f_S(0).
    fx0 <- fx[1]
    fx <- fx[-1]

    ## Distributions are expressed as a member of the (a, b, 0) or (a,
    ## b, 1) families of distributions. Assign parameters 'a' and 'b'
    ## depending of the chosen distribution and compute f_S(0) in
    ## every case, and p1 if p0 is specified in argument.
    dist <- model.freq$dist
    par <- model.freq$par
    if (dist == "geom")
    {
        dist <- "nbinom"
        par$size <- 1
    }

    if (dist == "pois")
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
    else if (dist == "nbinom")
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
    else if (dist == "binom")
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

    ## The recursion formula is slightly different for the (a, b, 0)
    ## and (a, b, 1) cases. We do the split here to avoid repeatedly
    ## testing in which case we're in.
    ##
    ## Vector 'fs' will hold the probabilities and will be expanded as
    ## needed. We are not supposed to do that in S, but assigning a
    ## longer than needed vector of NAs proved cumbersome and slower.
    fs <- fs0
    cumul <- sum(fs)

    ## (a, b, 0) case
    if (is.null(p0))
    {
        ## See in the (a, b, 1) case why this is defined here.
        r <- length(fx)

        repeat
        {
            if (echo)  print(tail(cumul, 1))

            x <- length(fs)
            m <- min(x, r)
            fs <- c(fs, sum((a + b * 1:m / x) * head(fx, m) * rev(tail(fs, m)))/(1 - a * fx0))
            if (TOL < (cumul <- cumul + tail(fs, 1)))
                break
        }
    }
    ## (a, b, 1) case
    else
    {
        ## Line below is a hack to reproduce the fact that the
        ## distribution of claim amounts is 0 past its maximum
        ## value. Only needed in the (a, b, 1) case for the additional
        ## term in the recursion formula.
        fx <- c(fx, 0)
        r <- length(fx)
        const <- p1 - (a + b) * p0

        repeat
        {
            if (echo) print(tail(cumul, 1))

            x <- length(fs)
            m <- min(x, r)
            fs <- c(fs, (const * fx[m] + sum((a + b * 1:m / x) * head(fx, m) * rev(tail(fs, m))))/(1 - a * fx0))
            if (TOL < (cumul <- cumul + tail(fs, 1)))
                break
        }
    }

    FUN <- stepfun((0:(length(fs) - 1)) * x.scale, c(0, cumsum(fs)))
    class(FUN) <- c("aggregateDist", "ecdf", class(FUN))
    assign("fs", fs, env = environment(FUN))
    assign("Call", Call, env = environment(FUN))
    assign("x.scale", x.scale, env = environment(FUN))
    comment(FUN) <- "Recursive method approximation"
    FUN
}
