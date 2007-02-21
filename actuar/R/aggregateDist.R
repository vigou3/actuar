### ===== actuar: an R package for Actuarial Science =====
###
### Use one of five methods to compute the aggregate claim amount
### distribution of a portfolio over a period given a frequency and a
### severity model or the true moments of the distribution.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Louis-Philippe Pouliot

aggregateDist <-
    function(method = c("recursive", "convolution", "normal", "npower", "simulation"),
             model.freq = NULL, model.sev = NULL, p0 = NULL, x.scale = 1,
             moments, nb.simul, ..., TOL = 1e-06, echo = FALSE)
{
    Call <- match.call()

    ## The method used essentially tells which function should be
    ## called for the calculation of the aggregate claims
    ## distribution.
    method <- match.arg(method)

    if (method == "normal")
    {
        ## An error message is issued if the number of moments listed
        ## is not appropriate for the method. However it is the user's
        ## responsability to list the moments in the correct order
        ## since the vector is not required to be named.
        if (missing(moments) || length(moments) < 2)
            stop("'moments' must supply the mean and variance of the distribution")
        return(normal(moments[1], moments[2]))
    }

    if (method == "npower")
    {
        if (missing(moments) || length(moments) < 3)
            stop("'moments' must supply the mean, variance and skewness of the distribution")
        return(npower(moments[1], moments[2], moments[3]))
    }

    if (method == "simulation")
    {
        if (missing(nb.simul))
            stop("'nb.simul' must supply the number of simulations")
        if (is.null(names(model.freq)) && is.null(names(model.sev)))
            stop("expressions in 'model.freq' and 'model.sev' must be named")
        return(simS(nb.simul, model.freq = model.freq, model.sev = model.sev))
    }

    ## "recursive" and "convolution" cases. Both require a discrete
    ## distribution of claim amounts, that is a vector of
    ## probabilities in argument 'model.sev'.
    if (!is.numeric(model.sev))
        stop("'model.sev' must be a vector of probabilities")

    ## Recursive method uses a model for the frequency distribution.
    if (method == "recursive")
    {
        if (is.null(model.freq) || !is.character(model.freq))
            stop("frequency distribution must be supplied as a character string")
        dist <- match.arg(tolower(model.freq),
                          c("poisson", "geometric", "negative binomial",
                            "binomial", "logarithmic"))
        return(panjer(fx = model.sev, dist = dist, p0 = p0,
                      x.scale = x.scale, ..., echo = echo, TOL = TOL))
    }

    ## Convolution method requires a vector of probabilites in
    ## argument 'model.freq'.
    if (method == "convolution")
    {
        if (!is.numeric(model.freq))
            stop("'model.freq' must be a vector of probabilities")
        return(exact(fx = model.sev, pn = model.freq, x.scale = x.scale))
    }

    stop("internal error")
}

print.aggregateDist <- function(x, ...)
{
    cat("\nAggregate Claim Amount Distribution\n")
    cat("  ", label <- comment(x), "\n\n")

    cat("Call:\n")
    print(get("Call", envir = environment(x)))
    cat("\n")

    if (label %in% c("Exact calculation (convolutions)",
                     "Recursive method approximation",
                     "Approximation by simulation"))
    {
        n <- length(get("x", environment(x)))
        cat("Data:  (", n, "obs. )\n")
        numform <- function(x) paste(formatC(x, dig = 4, width = 5), collapse = ", ")
        i1 <- 1:min(3, n)
        i2 <- if (n >= 4)
            max(4, n - 1):n
        else integer(0)
        xx <- eval(expression(x), env = environment(x))
        cat(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3)
        ", ", if (n > 5)
        " ..., ", numform(xx[i2]), "\n", sep = "")
        cat("\n")
    }
    if (label %in% c("Normal approximation",
                     "Normal Power approximation"))
        cat(attr(x, "source"), "\n")
    print(environment(x))
    cat("Class attribute:\n")
    print(attr(x, "class"))
}
