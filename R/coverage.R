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

coverage2 <- function(pdf, cdf, deductible = 0, franchise = FALSE,
                     limit = Inf, coinsurance = 1, inflation = 0,
                     per.loss = FALSE)
{
    Call <- match.call()

    ## First determine if the cdf is needed or not. It is needed when
    ## there is a deductible or a limit and, of course, if the output
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

    ## We will build the body of the output function as an expression,
    ## piece by piece as we go along. We start here with the first
    ## couple of lines that never change.
    e <- expression(
        args <- as.list(match.call())[-c(1L, 2L)],
        nargs <- names(args))

    ## Prepare the cdf object for cases needing the cdf.
    if (needs.cdf)
    {
        ## Argument 'cdf' can be a character string giving the
        ## function name or a straight function object. In the latter
        ## case, 'cdf' will contain function definitions and the
        ## output function of coverage() will contain multiple
        ## function definitions. Not nice. Instead, always treat 'cdf'
        ## (and 'pdf', below) as character strings.
        cdf <- as.character(Call$cdf)

        ## Initialization in the output function.
        e <- c(e,
               substitute(F <- fun, list(fun = as.name(cdf))),
               quote(formals(F)[nargs] <- args))

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

        ## Definitions of 1 - F(d) and 1 - F(u), using 'lower.tail =
        ## FALSE' if available in 'cdf'.
        if (has.lower)
        {
            Sd <- substitute(F(a, lower.tail = FALSE), list(a = d))
            Su <- substitute(F(a, lower.tail = FALSE), list(a = u))
        }
        else
        {
            Sd <- substitute(1 - F(a), list(a = d))
            Su <- substitute(1 - F(a), list(a = u))
        }
    }

    ## Repeat same steps as above for case needing the pdf: output
    ## function is a pdf.
    if (!is.cdf)
    {
        pdf <- as.character(Call$pdf)   # same as 'cdf' above
        e <- c(e,
               substitute(f <- fun, list(fun = as.name(pdf))),
               quote(formals(f)[nargs] <- args))
        formalsPDF <- formals(pdf)      # arguments as list
        argsPDF <- setdiff(names(formalsPDF), "log") # drop argument 'log'
        argsFUN <- formalsPDF[argsPDF]  # arguments of output function
        x <- as.name(argsPDF[1])        # symbol
    }

    ## Initialization of the results vector in the output function
    ## with 0s.
    e <- c(e,
           substitute(res <- numeric(length(x)), list(x = x)))

    ## Build the value at which the underlying pdf/cdf will be called
    ## for non special case values of 'x'.
    x.mod <- as.call(c(as.name("["), x, as.name("w")))
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
            quote(x == 0)
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

    ## Definition for the first branch. Computation to make only if
    ## there is a decutible with the payment per loss random variable.
    ## For all other cases, the value in the first branch is 0.
    if (per.loss & deductible)
        e <- c(e,
               substitute(w <- which(cond1), list(cond1 = cond1)),
               substitute(res[w] <- F(a), list(a = d)))

    ## Function definitions for the second and third branches. The
    ## 'is.cdf = TRUE' and 'is.cdf = FALSE' cases must be treated
    ## separately.
    if (is.cdf)
    {
        cond3 <- substitute(x >= b, list(x = x, b = bound2))
        f2 <- substitute(F(a), list(a = x.mod))
        f3 <- 1
        if (!per.loss & deductible)
            f2 <- substitute((f - F(a))/S,
                             list(f = f2, S = Sd, a = d))
    }
    else
    {
        cond3 <- substitute(x == b, list(x = x, b = bound2))
        f2 <- substitute(f(a), list(a = x.mod))
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

    ## Definitions for the second and third branch in the output
    ## function.
    e <- c(e,
           substitute(w <- which(cond2), list(cond2 = cond2)),
           substitute(res[w] <- f2, list(f2 = f2)),
           substitute(w <- which(cond3), list(cond3 = cond3)),
           substitute(res[w] <- f3, list(f3 = f3)),
           quote(res))

    ## Wrap up the output function.
    body(FUN) <- as.call(c(as.name("{"), e)) # taken from help(body)
    formals(FUN) <- argsFUN             # set arguments
    environment(FUN) <- new.env()       # new, empty environment
    FUN
}
