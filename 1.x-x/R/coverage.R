### ===== actuar: An R Package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution
### function for data with deductible, limit, coinsurance and
### inflation.
###
### See Chapter 5 of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

coverage <- function(pdf, cdf, deductible = 0, franchise = FALSE,
                     limit = Inf, coinsurance = 1, inflation = 0,
                     per.loss = FALSE)
{
    Call <- match.call()

    ## First determine if the cdf is needed or not. It is when there
    ## is a deductible or a limit and, of course, if the output
    ## function should compute the cdf.
    is.cdf <- missing(pdf) || is.null(pdf) # return cdf?
    needs.cdf <- any(deductible > 0, limit < Inf, is.cdf) # cdf needed?

    ## Sanity check of arguments
    if (any(deductible < 0, limit < 0, coinsurance < 0, inflation < 0))
        stop("coverage modifications must be positive")
    if (limit <= deductible)
      stop("deductible must be smaller than the limit")
    if (coinsurance > 1)
        stop("coinsurance must be between 0 and 1")
    if (missing(cdf) & needs.cdf)
        stop("'cdf' must be supplied")

    ## Quantites often used
    r <- 1 + inflation
    d <- deductible/r
    u <- limit/r

    ## Prepare the cdf object for cases needing the cdf: output
    ## function is a cdf or there is a deductible or a limit.
    if (needs.cdf)
    {
        ## Arguments 'cdf' can be a character string giving the
        ## function name or a straight function object. In the latter
        ## case, 'cdf' will contain function definitions and the
        ## output function of coverage() will contain multiple
        ## function definitions. Not nice. Instead, always treat 'cdf'
        ## (and 'pdf', below) as character strings.
        cdf <- as.character(Call$cdf)

        ## Get argument list to build function calls and, eventually,
        ## to specify arguments of the output function.
        formalsCDF <- formals(cdf)          # arguments as list
        argsCDF <- names(formalsCDF)        # arguments names as strings

        ## Remember if argument 'lower.tail' is available, so we can
        ## use it later. Then, drop unsupported arguments 'lower.tail'
        ## and 'log.p'.
        has.lower <- "lower.tail" %in% argsCDF
        argsCDF <- setdiff(argsCDF, c("lower.tail", "log.p"))

        ## If output function is a cdf:
        ##
        ## 1. Set its arguments to those of 'cdf'.
        ## 2. Set the symbol representing the variable in function
        ##    calls. Should be the first argument of the output
        ##    function.
        ## 3. Drop the first argument of 'cdf' since it is no longer
        ##    used after this block.
        if (is.cdf)
        {
            argsFUN <- formalsCDF[argsCDF]  # arguments of output function
            x <- as.name(argsCDF[1])        # symbol
        }

        ##  Prepare argument list for use in do.call
        argsCDF <- sapply(argsCDF[-1], as.name) # for use in do.call()

        ## Definitions of 1 - F(d) and 1 - F(u), using 'lower.tail =
        ## FALSE' if available in 'cdf'.
        if (has.lower)
        {
            Sd <- substitute(do.call(F, a),
                             list(F = cdf, a = c(d, argsCDF, lower.tail = FALSE)))
            Su <- substitute(do.call(F, a),
                             list(F = cdf, a = c(u, argsCDF, lower.tail = FALSE)))
        }
        else
        {
            Sd <- substitute(1 - do.call(F, a),
                             list(F = cdf, a = c(d, argsCDF)))
            Su <- substitute(1 - do.call(F, a),
                             list(F = cdf, a = c(u, argsCDF)))
        }
    }

    ## Repeat same steps as above for case needing the pdf: output
    ## function is a pdf.
    if (!is.cdf)
    {
        pdf <- as.character(Call$pdf)   # same as 'cdf' above
        formalsPDF <- formals(pdf)      # arguments as list
        argsPDF <- setdiff(names(formalsPDF), "log") # drop argument 'log'
        argsFUN <- formalsPDF[argsPDF]  # arguments of output function
        x <- as.name(argsPDF[1])        # symbol
        argsPDF <- sapply(argsPDF[-1], as.name) # for use in do.call()
    }

    ## Build the value at which the underlying pdf/cdf will be called
    ## for non special case values of 'x'.
    x.mod <- x
    if (coinsurance < 1)
        x.mod <- substitute(x/alpha, list(x = x.mod, alpha = coinsurance))
    if (deductible & !franchise)
        x.mod <- substitute(x + d, list(x = x.mod, d = deductible))
    if (inflation)
        x.mod <- substitute((x)/r, list(x = x.mod, r = r))

    ## Each pdf/cdf is defined in three branches. Define the
    ## boundaries and conditions for the first two branches.
    if (franchise)
    {
        bound1 <- coinsurance * deductible
        bound2 <- coinsurance * limit
        cond1 <- if (is.cdf)
            substitute(0 <= x & x <= b1, list(x = x, b1 = bound1))
        else
            substitute(x == 0, list(x = x))
        cond2 <- substitute(b1 < x & x < b2,
                            list(x = x, b1 = bound1, b2 = bound2))
    }
    else
    {
        bound1 <- 0
        bound2 <- coinsurance * (limit - deductible)
        cond1 <- substitute(x == 0, list(x = x))
        cond2 <- substitute(0 < x & x < b, list(x = x, b = bound2))
    }

    ## Function definition for the first branch.
    f1 <- if (per.loss & deductible)
        substitute(do.call(F, a), list(F = cdf, a = c(d, argsCDF)))
    else 0

    ## Function definitions for the second and third branches. The
    ## 'is.cdf = TRUE' and 'is.cdf = FALSE' cases must be treated
    ## separately.
    if (is.cdf)
    {
        cond3 <- substitute(x >= b, list(x = x, b = bound2))
        f2 <- substitute(do.call(F, a),
                         list(F = cdf, a = c(x.mod, argsCDF)))
        f3 <- 1
        if (!per.loss & deductible)
            f2 <- substitute((f - do.call(F, d))/S,
                             list(f = f2, F = cdf, S = Sd, d = c(d, argsCDF)))
    }
    else
    {
        cond3 <- substitute(x == b, list(x = x, b = bound2))
        f2 <- substitute(do.call(f, a),
                         list(f = pdf, a = c(x.mod, argsPDF)))
        f3 <- if (is.finite(limit)) Su else 0
        if (!per.loss & deductible)
        {
            f2 <- substitute(f/S, list(f = f2, S = Sd))
            if (is.finite(limit))
                f3 <- substitute(f/S, list(f = f3, S = Sd))
        }
        if (inflation | coinsurance < 1)
            f2 <- substitute(f/k, list(f = f2, k = coinsurance * r))
    }

    ## Output function
    eval(substitute(FUN <- function()
                    ifelse(cond1, f1,
                           ifelse(cond2, f2,
                                  ifelse(cond3, f3, 0))),
                    list(cond1 = cond1, cond2 = cond2, cond3 = cond3,
                         f1 = f1, f2 = f2, f3 = f3)))
    formals(FUN) <- argsFUN             # set arguments
    environment(FUN) <- new.env()       # new, empty environment
    FUN
}
