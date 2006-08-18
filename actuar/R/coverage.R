### ===== actuar: an R package for Actuarial Science =====
###
### Create modified density and modified cumulative distribution function
### for data with deductible d and limit u.
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

coverage <- function(dist, param, d, u, density = TRUE, per.payment = TRUE, ...)
{
    f1 <- match.fun(paste("d", dist$dist, sep = ""))
    formals(f1)[dist$par] <- param
    F1 <- match.fun(paste("p", dist$dist, sep = ""))
    formals(F1)[dist$par] <- param
    Call <- match.call()
    if (any(u < d))
      stop("deductible must be less than limit")
    
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

