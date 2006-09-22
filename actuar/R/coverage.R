### ===== actuar: an R package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution
### function for data with deductible, limit, coinsurance and
### inflation.
###
### See Chapter 5 of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

coverage <- function(dist, deductible = 0, franchise = FALSE,
                     limit = Inf, coinsurance = 1, inflation = 0,
                     per.loss = FALSE, cdf = FALSE)
{
    ## Sanity check of arguments
    if (any(deductible < 0, limit < 0, coinsurance < 0, inflation < 0))
        stop("coverage modifications must be positive")
    if (limit <= deductible)
      stop("deductible must be smaller than the limit")
    if (coinsurance > 1)
        stop("coinsurance must be between 0 and 1")

    ## Function 'pdist'() is always needed. Get its argument list to
    ## build function calls and, eventually, specify arguments of the
    ## output function.
    F <- paste("p", dist, sep = "")     # 'pdist'()
    formalsF <- formals(F)              # arguments as list
    argsF <- names(formalsF)            # arguments names as strings

    ## Remember if argument 'lower.tail' is available, so we can use
    ## it later. Then, drop unsupported arguments 'lower.tail' and
    ## 'log.p'.
    has.lower <- if ("lower.tail" %in% argsF) TRUE else FALSE
    argsF <- setdiff(argsF, c("lower.tail", "log.p"))

    ## 1. Set arguments of the output function. Should be those of
    ##    'pdist'() or 'ddist'() depending if 'cdf' is 'TRUE' or
    ##    'FALSE', respectively.
    ##
    ## 2. Set the symbol representing the variable in function
    ##    calls. Should be the first argument of the output
    ##    function.
    ##
    ## 3. Drop unsupported argument 'log', if present, in 'ddist'().
    ##
    ## 4. Drop the first argument of 'pdist'() and 'ddist'() which are
    ##    no longer used after this block and prepare argument list
    ##    for use in do.call().
    if (cdf)
    {
        argsFUN <- formalsF[argsF]      # arguments of output function
        x <- as.name(argsF[1])          # symbol
    }
    else
    {
        f <- paste("d", dist, sep = "") # 'ddist'()
        formalsf <- formals(f)          # arguments as list
        argsf <- setdiff(names(formalsf), "log") # drop argument 'log'
        argsFUN <- formalsf[argsf]      # arguments of output function
        x <- as.name(argsf[1])          # symbol
        argsf <- sapply(argsf[-1], as.name) # for use in do.call()
    }
    argsF <- sapply(argsF[-1], as.name) # for use in do.call()

    ## Quantites often used
    r <- 1 + inflation
    d <- deductible/r
    u <- limit/r

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
        cond1 <- if (cdf)
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

    ## Definitions of 1 - F(d) and 1 - F(u), using 'lower.tail =
    ## FALSE' if available in 'pdist'().
    if (has.lower)
    {
        Sd <- substitute(do.call(F, a),
                         list(F = F, a = c(d, argsF, lower.tail = FALSE)))
        Su <- substitute(do.call(F, a),
                         list(F = F, a = c(u, argsF, lower.tail = FALSE)))
    }
    else
    {
        Sd <- substitute(1 - do.call(F, a),
                         list(F = F, a = c(d, argsF)))
        Su <- substitute(1 - do.call(F, a),
                         list(F = F, a = c(u, argsF)))
    }

    ## Function definition for the first branch.
    f1 <- if (per.loss & deductible)
              substitute(do.call(F, a), list(F = F, a = c(d, argsF)))
          else 0

    ## Function definitions for the second and third branches. The
    ## 'cdf = TRUE' and 'CDF = FALSE' must be treated separately.
    if (cdf)
    {
        cond3 <- substitute(x >= b, list(x = x, b = bound2))
        f2 <- substitute(do.call(F, a),
                         list(F = F, a = c(x.mod, argsF)))
        f3 <- 1
        if (!per.loss & deductible)
            f2 <- substitute((f - do.call(F, d))/S,
                             list(f = f2, F = F, S = Sd, d = c(d, argsF)))
    }
    else
    {
        cond3 <- substitute(x == b, list(x = x, b = bound2))
        f2 <- substitute(do.call(f, a),
                         list(f = f, a = c(x.mod, argsf)))
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
