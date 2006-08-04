### ===== actuar: an R package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution function
### for data with deductible d and limit u.
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>


distrmod <- function(dist, param, d, u, density = TRUE, per.payment = TRUE, ...)
{
    f1 <- match.fun(paste("d", dist$dist, sep = ""))
    formals(f1)[dist$par] <- param
    F1 <- match.fun(paste("p", dist$dist, sep = ""))
    formals(F1)[dist$par] <- param
    Call <- match.call()
    if (u < d)
      stop("deductible must be less than limit.")
    
    ## Create modified density
    if (density)
    {
        if (per.payment)
        {
            FUN <- function(x, ...)
                ifelse (x >= 0 & x < (u - d), f1(x + d) / (1 - F1(d)), ifelse (x == (u - d), (1 - F1(u)) / (1 - F1(d)), 0))
        }
        else
        {
            FUN <- function(x, ...)
                ifelse (x == 0, F1(d), ifelse (x > 0 & x < (u - d), f1(x + d), ifelse (x == (u - d), 1 - F1(u), 0)))
        }
    }
    ## Create modified cumulative distribution function
    else
    {
        if (per.payment)
        {
            FUN <- function(x, ...)
                ifelse (x >= 0 & x < (u - d), (F1(x + d) - F1(d)) / (1 - F1(d)), 1)
        }
        else
        {
            FUN <- function(x, ...)
                ifelse (x == 0, F1(d), ifelse(x > 0 & x < (u - d), F1(x + d), 1))
        }
    }
    attr(FUN, "call") <- Call
    FUN
}

mdistrmod <- function(dist, param, d, u, a = 1, r = 0, ...)
{
    if (u < d)
       stop("deductible must be less than limit.")
    if (a < 0 || a > 1)
      stop("coinsurance must be between 0 and 1.")
    lev1 <- match.fun(paste("lev", dist$dist, sep = ""))
    formals(lev1)[dist$par] <- param
    lev2 <- match.fun(paste("lev", dist$dist, sep = ""))
    formals(lev2)[dist$par] <- param
    formals(lev2)["order"] <- 2
    F1 <- match.fun(paste("p", dist$dist, sep = ""))
    formals(F1)[dist$par] <- param
    Call <- match.call()

    Eloss <- a * (1 + r) * (lev1(u / (1 + r)) - lev1(d / (1 + r)))
    Epmt <- Eloss / (1 - F1(d / (1 + r)))
    E2loss <- a ^ 2 * (1 + r) ^ 2 * (lev2(u / (1 + r)) - lev2(d / (1 + r)) - 2 * (d / (1 + r)) * lev1(u / (1 + r)) + 2 * (d / (1 + r)) * lev1(d / (1 + r)))
    E2pmt <- E2loss / (1 - F1(d / (1 + r)))
    Vloss <- E2loss - Eloss ^ 2
    Vpmt <- E2pmt - Epmt ^ 2
    res = list(Eloss = Eloss, Epmt = Epmt, Varloss = Vloss, Varpmt = Vpmt)
    class(res) <- c("mdistrmod", class(res))
    attr(res, "call") <- sys.call()
    res
}

print.mdistrmod <- function(x, ...)
{
    cat("Expected-value of the per-payment variable is :", x$Eloss, "\n",
    "Variance of the per-payment variable is :", x$Varloss, "\n",
    "Expected-value of the per-loss variable is :", x$Epmt, "\n",
    "Variance of the per-loss variable is :", x$Varpmt, "\n")
}
