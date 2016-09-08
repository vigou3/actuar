### ===== actuar: An R Package for Actuarial Science =====
###
### Panjer recursion formula to compute the approximate aggregate
### claim amount distribution of a portfolio over a period.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sebastien Auclair, Louis-Philippe Pouliot and Tommy Ouellet

panjer <- function(fx, dist, p0 = NULL, x.scale = 1, ...,
                   convolve = 0, tol = sqrt(.Machine$double.eps),
                   maxit = 500, echo = FALSE)
{
    ## Express 'tol' as a value close to 1. If needed, modify the
    ## accuracy level so that the user specified level is attained
    ## *after* the additional convolutions (without getting too high).
    tol <- if (convolve > 0)
               min((0.5 - tol + 0.5)^(0.5 ^ convolve),
                   0.5 - sqrt(.Machine$double.eps) + 0.5)
           else
               0.5 - tol + 0.5

    ## Check if p0 is a valid probability.
    if (!is.null(p0))
    {
        if (length(p0) > 1L)
        {
            p0 <- p0[1L]
            warning("'p0' has many elements: only the first used")
        }
        if ((p0 < 0) || (p0 > 1))
            stop("'p0' must be a valid probability (between 0 and 1)")
    }

    ## Treat trivial case where 'p0 == 1' and hence F_S(0) = 1.
    if (identical(p0, 1))
    {
        FUN <- approxfun(0, 1, method = "constant",
                         yleft = 0, yright = 1, f = 0)
        class(FUN) <- c("ecdf", "stepfun", class(FUN))
        assign("fs", 1, envir = environment(FUN))
        assign("x.scale", x.scale, envir = environment(FUN))
        return(FUN)
    }

    ## The call to .External below requires 'p1' to be initialized.
    p1 <- 0

    ## Argument '...' should contain the values of the parameters of
    ## 'dist'.
    par <- list(...)

    ## Distributions are expressed as a member of the (a, b, 0) or (a,
    ## b, 1) families of distributions. Assign parameters 'a' and 'b'
    ## depending of the chosen distribution and compute f_S(0) in
    ## every case, and p1 if p0 is specified in argument.
    ##
    ## At this point, either p0 is NULL or 0 <= p0 < 1.
    if (startsWith(dist, "zero truncated"))
    {
        if ( !(is.null(p0) || identical(p0, 0)) )
            warning("value of 'p0' is ignored with a zero truncated distribution")
        dist = sub("zero truncated ", "", dist) # drop "zero truncated" prefix
        p0 <- 0
    }

    if (dist == "geometric")
    {
        dist <- "negative binomial"
        par$size <- 1
    }

    if (dist == "poisson")
    {
        if (!"lambda" %in% names(par))
            stop("value of 'lambda' missing")
        lambda <- par$lambda
        a <- 0
        b <- lambda
        if (is.null(p0))
            fs0 <- exp(lambda * (fx[1L] - 1))
        else  # 0 <= p0 < 1
        {
            fs0 <- p0 + (1 - p0) * pgfztpois(fx[1L], lambda)
            p1 <- (1 - p0) * dztpois(1, lambda)
        }
    }
    else if (dist == "negative binomial")
    {
        if (!all(c("prob", "size") %in% names(par)))
            stop("value of 'prob' or 'size' missing")
        p <- par$prob
        r <- par$size
        a <- 1 - p
        b <- (r - 1) * a
        if (is.null(p0))
            fs0 <- exp(-r * log1p(-a/p * (fx[1L] - 1)))
        else  # 0 <= p0 < 1
        {
            fs0 <- p0 + (1 - p0) * pgfztnbinom(fx[1L], r, p)
            p1 <- (1 - p0) * dztnbinom(1, r, p)
        }
    }
    else if (dist == "binomial")
    {
        if (!all(c("prob", "size") %in% names(par)))
            stop("value of 'prob' or 'size' missing")
        a <- - p/(1 - p)
        b <- -(n + 1) * a
        if (is.null(p0))
            fs0 <- exp(n * log1p(p * (fx[1L] - 1)))
        else  # 0 <= p0 < 1
        {
            fs0 <- p0 + (1 - p0) * pgfztbinom(fx[1L], n, p)
            p1 <- (1 - p0) * dztbinom(1, n, p)
        }
    }
    else if (dist == "logarithmic")
    {
        if (!"prob" %in% names(par))
            stop("value of 'prob' missing")
        a <- prob
        b <- -a
        if (is.null(p0) || identical(p0, 0)) # logarithmic always "zero truncated"
            fs0 <- pgflogarithmic(fx[1L], prob)
        else # 0 < p0 < 1; quite unusual "zero modified logarithmic"
        {
            fs0 <- p0 + (1 - p0) * pgflogarithmic(fx[1L], prob)
            p1 <- (1 - p0) * dlogarithmic(1, prob)
        }
    }
    else
        stop("frequency distribution not in the (a, b, 0) or (a, b, 1) families")

    ## If fs0 is equal to zero, the recursion will not start. There is
    ## no provision to automatically cope with this situation in the
    ## current version of this function. Just issue an error message
    ## and let the user do the work by hand.
    if (identical(fs0, 0))
        stop("Pr[S = 0] is numerically equal to 0; impossible to start the recursion")

    ## Recursive calculations in C.
    fs <- .External("actuar_do_panjer", p0, p1, fs0, fx, a, b, convolve, tol, maxit, echo)

    FUN <- approxfun((0:(length(fs) - 1)) * x.scale, pmin(cumsum(fs), 1),
                     method = "constant", yleft = 0, yright = 1, f = 0,
                     ties = "ordered")
    class(FUN) <- c("ecdf", "stepfun", class(FUN))
    assign("fs", fs, envir = environment(FUN))
    assign("x.scale", x.scale, envir = environment(FUN))
    FUN
}
