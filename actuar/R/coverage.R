### ===== actuar: an R package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution function
### for data with deductible d and limit u.
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

coverage <- function(dist, deductible = 0, franchise = FALSE,
                     limit = Inf, coinsurance = 1, inflation = 0,
                     per.loss = FALSE, cdf = FALSE)
{
    if (any(deductible < 0, limit < 0, coinsurance < 0, inflation < 0))
        stop("coverage modifications must be positive")
    if (limit <= deductible)
      stop("deductible must be smaller than the limit")

    if (!cdf)
        f <- paste("d", dist, sep = "")
    F <- paste("p", dist, sep = "")
    args <- sapply(names(formals(F)[-1]), as.name)

    r <- 1 + inflation
    d <- deductible/r
    u <- limit/r

    y <- as.name("x")

    if (coinsurance < 1)
        y <- substitute(x/alpha, list(x = y, alpha = coinsurance))
    if (deductible & !franchise)
        y <- substitute(x + d, list(x = y, d = deductible))
    if (inflation)
        y <- substitute((x)/r, list(x = y, r = r))

    if (franchise)
    {
        bound1 <- coinsurance * deductible
        bound2 <- coinsurance * limit
        cond1 <- substitute(0 <= x & x <= b1, list(b1 = bound1))
        cond2 <- substitute(b1 < x & x < b2, list(b1 = bound1, b2 = bound2))
    }
    else
    {
        bound1 <- 0
        bound2 <- coinsurance * (limit - deductible)
        cond1 <- substitute(x == 0)
        cond2 <- substitute(0 < x & x < b, list(b = bound2))
    }

    if (per.loss)
        f1 <- substitute(do.call(F, a), list(F = F, a = c(x = d, args)))
    else
        f1 <- 0

    if (cdf)
    {
        if (per.loss)
            f2 <- substitute(do.call(F, a), list(F = F, a = c(x = y, args)))
        else
            f2 <- substitute((do.call(F, a) - do.call(F, d))/(1 - do.call(F, d)),
                             list(F = F, a = c(x = y, args), d = c(x = d, args)))
        cond3 <- substitute(x >= b, list(b = bound2))
        f3 <- 1
    }
    else
    {
        f2 <- substitute(do.call(f, a), list(f = f, a = c(x = y, args)))
        if (!per.loss)
            f2 <- substitute(f/(1 - do.call(F, d)),
                             list(f = f2, F = F, d = c(x = d, args)))
        if (inflation | coinsurance < 1)
            f2 <- substitute(f/k,
                             list(f = f2, k = coinsurance * r))

        cond3 <- substitute(x == b, list(b = bound2))
        f3 <- substitute(1 - do.call(F, a), list(F = F, a = c(x = u, args)))
        if (!per.loss)
            f3 <- substitute(f/(1 - do.call(F, d)),
                             list(f = f3, F = F, d = c(x = d, args)))
    }

    eval(substitute(FUN <- function()
               ifelse(cond1, f1,
                      ifelse(cond2, f2,
                             ifelse(cond3, f3, 0))),
               list(cond1 = cond1, cond2 = cond2, cond3 = cond3,
                    f1 = f1, f2 = f2, f3 = f3)))
    formals(FUN) <- formals(pdf)
    FUN
}


### TESTS

## Franchise ordinaire seulement, par sinitre
coverage("gamma", deductible = 1000, per.loss = TRUE, cdf = TRUE)
coverage("gamma", deductible = 1000, per.loss = TRUE, cdf = FALSE)

## Tout le toutim, franchise ordinaire, par sinistre
coverage("gamma", deductible = 1, limit = 10, coinsurance = 0.9,
         inflation = 0.05, per.loss = TRUE, cdf = TRUE)
