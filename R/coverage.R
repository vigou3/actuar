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

    ## Quantites often used.
    r <- 1 + inflation
    d <- deductible/r
    u <- limit/r

    ## The modified cdf or pdf are usually defined in branches. To
    ## avoid using nested ifelse(), we will rather rely on sets of
    ## expressions to make the required calculations for each branch
    ## separately. This is actually much faster.
    ##
    ## The output function will have varying number of said
    ## expressions depending on the case that is dealt with. We will
    ## build the body of the output function piece by piece as we go
    ## along. We start here with the first few lines that never
    ## change.
    e <- expression(Call <- match.call())

    ## One main discriminating factor is whether the cdf is needed for
    ## the output function of not.
    if (needs.cdf)
    {
        ## Get argument list of 'cdf' to transfert them to the output
        ## function.
        argv <- formals(cdf)            # arguments as list
        argn <- names(argv)             # arguments names as strings

        ## Remember if argument 'lower.tail' is available, so we can
        ## use it later. Then, drop unsupported arguments 'lower.tail'
        ## and 'log.p'.
        has.lower <- "lower.tail" %in% argn
        argn <- setdiff(argn, c("lower.tail", "log.p"))

        ## If output function is a cdf:
        ##
        ## 1. Set its arguments to those of 'cdf'.
        ## 2. Set the symbol representing the variable in function
        ##    calls. Should be the first argument of the output
        ##    function.
        ## Initialization of the output function.
        if (is.cdf)
        {
            ## The output function will have the cdf formal arguments.
            ## Object 'x' holds the symbol of the first argument.
            argsFUN <- argv[argn]       # arguments of output function
            x <- as.name(argn[1])       # symbol

            ## A function to compute F(y) is needed; one to compute
            ## F(d) only if there is a deductible; one to also compute
            ## 1 - F(d) for the per payment cases. Never need to
            ## compute F(u).
            e <- c(e,
                   quote(F <- Call),
                   substitute(F[[1L]] <- as.name(fun),
                              list(fun = as.character(Call$cdf))))
            if (deductible)
            {
                e <- c(e,
                       Fd <- F,
                       substitute(Fd[[2L]] <- a, list(a = d)))
                if (!per.loss & has.lower)
                    e <- c(e,
                           Sd <- Fd,
                           quote(Sd$lower.tail <- FALSE))
            }
        }
        else
        {
            ## No function to compute F(y) needed. However, need to
            ## compute F(d) for the per loss cases and 1 - F(d) for
            ## the per payment cases; compute 1 - F(u) when there is a
            ## limit.
            if (deductible)
            {
                Fd.name <- as.name(if (!per.loss & has.lower) "Sd" else "Fd")
                e <- c(e,
                       substitute(G <- Call, list(G = Fd.name)),
                       substitute(G[[1L]] <- as.name(fun),
                                  list(G = Fd.name, fun = as.character(Call$cdf))),
                       substitute(G[[2L]] <- a, list(G = Fd.name, a = d)),
                       if (!per.loss & has.lower) quote(Sd$lower.tail <- FALSE))
                if (limit < Inf)
                {
                    Fu.name <- as.name(if (has.lower) "Su" else "Fu")
                    e <- c(e,
                           substitute(H <- G, list(H = Fu.name, G = Fd.name)),
                           substitute(Su[[2L]] <- a, list(a = u)),
                           if (per.loss & has.lower) quote(Su$lower.tail <- FALSE))
                }
            }
            else
            {
                ## Since needs.cdf == TRUE, then this case necessarily
                ## has limit < Inf and we need to compute F(u).
                e <- c(e,
                       Su <- Call,
                       substitute(Su[[1L]] <- as.name(fun),
                                  list(fun = as.character(Call$cdf))),
                       substitute(Su[[2L]] <- a, list(a = u)),
                       if (has.lower) quote(Su$lower.tail <- FALSE))
            }
        }

        ## Fd <- Call,
        ## substitute(Fd[[1L]] <- fun, list(fun = as.name(Call$cdf))),
        ## substitute(names(Fd)[[2L]] <- q, list(q = argn[1L])),
        ## if (has.lower) Fd$lower.tail <- FALSE,
        ## substitute(Fd[[2L]] <- x, list(x = d)),
        ## Fu <- Fd,
        ## Fu[[2L] <- x, list(x = u))

        ## We will need to compute 1 - F(d) and 1 - F(u) in future
        ## expressions. We definite the expressions to do these
        ## computations here; they will differ depending if
        ## 'lower.tail = FALSE' is available or not in 'cdf'.
        ## if (has.lower)
        ## {
        ##     Fd$lower.tail <- FALSE
        ##     Sd <- substitute(Fd(a), list(a = d))
        ##     Su <- substitute(Fu(a), list(a = u))
        ## }
        ## else
        ## {
        ##     Sd <- substitute(1 - Fd(a), list(a = d))
        ##     Su <- substitute(1 - Fu(a), list(a = u))
        ## }
    }

    ## Repeat same steps as above for case needing the pdf. The output
    ## function is a pdf and in this case the arguments of the output
    ## function are those of 'pdf'.
    if (!is.cdf)
    {
        argv <- formals(pdf)                # arguments as list
        argn <- setdiff(names(argv), "log") # drop argument 'log'
        argsFUN <- argv[argn]           # arguments of output function
        x <- as.name(argn[1])           # symbol
        print(mode(e))
        e <- c(e,
               quote(f <- Call),
               substitute(f[[1L]] <- as.name(fun),
                          list(fun = as.character(Call$pdf))))
    }

    ## Build the value at which the underlying pdf/cdf will be called
    ## for non special case values of 'x'. We need to index 'x' to
    ## only compute for the correct values of a given branch.
    x.mod <- as.call(c(as.name("["), x, as.name("w")))
    if (coinsurance < 1)
        x.mod <- substitute(x/alpha, list(x = x.mod, alpha = coinsurance))
    if (deductible & !franchise)
        x.mod <- substitute(x + d, list(x = x.mod, d = deductible))
    if (inflation)
        x.mod <- substitute((x)/r, list(x = x.mod, r = r))

    ## Each pdf/cdf is defined in three branches. Define the
    ## boundaries and conditions for the first two branches. Those for
    ## the third branch are defined further down.
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

    ## Initialization of the results vector in the output function
    ## with 0s.
    e <- c(e,
           substitute(res <- numeric(length(x)), list(x = x)))

    ## Definition of the output function for the first branch. There
    ## is a computation to make only if there is a deductible with the
    ## payment per loss random variable. For all other cases, the
    ## value in the first branch is 0 and we rely on the
    ## initialization with numeric() done at the previous step.
    if (per.loss & deductible)
        e <- c(e,
               substitute(res[which(cond1)] <- eval.parent(Fd),
                          list(cond1 = cond1)))

    ## Definition of the output function for the second and third
    ## branches. The 'is.cdf = TRUE' and 'is.cdf = FALSE' cases must
    ## be treated separately.
    if (is.cdf)
    {
        cond3 <- substitute(x >= b, list(x = x, b = bound2))
        f2 <- quote(eval.parent(F))
        if (!per.loss & deductible)
            f2 <- substitute((f - g)/h,
                             list(f = f2,
                                  g = if (has.lower) quote(eval.parent(Fd)) else quote((q <- eval.parent(Fd))),
                                  h = if (has.lower) quote(eval.parent(Sd)) else quote((1 - q))))
        e <- c(e,
               substitute(w <- which(cond), list(cond = cond2)),
               substitute(F[[2L]] <- x, list(x = x.mod)),
               substitute(res[w] <- f, list(f = f2)),
               substitute(res[cond] <- 1, list(cond = cond3)))
    }
    else
    {
        cond3 <- substitute(x == b, list(x = x, b = bound2))
        Sd.quote <- if (has.lower) quote(eval.parent(Sd)) else quote((1 - eval.parent(Fd)))
        f2 <- quote(eval.parent(f))
        f3 <- if (has.lower) quote(eval.parent(Su)) else quote((1 - eval.parent(Fu)))
        if (!per.loss & deductible)
        {
            f2 <- substitute(f/S, list(f = f2, S = Sd.quote))
            if (limit < Inf)
                f3 <- substitute(f/S, list(f = f3, S = Sd.quote))
        }
        if (inflation | coinsurance < 1)
            f2 <- substitute(f/k, list(f = f2, k = coinsurance * r))
        e <- c(e,
               substitute(w <- which(cond), list(cond = cond2)),
               substitute(f[[2L]] <- x, list(x = x.mod)),
               substitute(res[w] <- f, list(f = f2)),
               substitute(res[cond] <- f, list(cond = cond3, f = f3)))
    }

    ## Definitions for the second and third branch in the output
    ## function.
    ## e <- c(e,
    ##        substitute(wx <- which(cond2), list(cond2 = cond2)),
    ##        substitute(res[wx] <- f2, list(f2 = f2)),
    ##        substitute(res[which(cond3)] <- f3, list(cond3 = cond3, f3 = f3)),
    ##        quote(res))

    ## Wrap up the output function.
    FUN <- function() {}
    body(FUN) <- as.call(c(as.name("{"), e)) # taken from help(body)
    formals(FUN) <- argsFUN             # set arguments
    environment(FUN) <- new.env()       # new, empty environment
    FUN
}
