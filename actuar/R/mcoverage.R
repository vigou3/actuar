### ===== actuar: an R package for Actuarial Science =====
###
### Calculate first and second moments of modified distribution 
### for data with deductible d and limit u.
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>


mcoverage <- function(dist, param, d, u, a = 1, r = 0, ...)
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

print.mcoverage <- function(x, ...)
{
    cat("Expected-value of the per-payment variable is :", x$Eloss, "\n",
    "Variance of the per-payment variable is :", x$Varloss, "\n",
    "Expected-value of the per-loss variable is :", x$Epmt, "\n",
    "Variance of the per-loss variable is :", x$Varpmt, "\n")
}
