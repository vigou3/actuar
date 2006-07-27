### ===== actuar: an R package for Actuarial Science =====
###
### Sample empirical limited value functions for individual and
### grouped data.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

elev <- function(x, ...)
    UseMethod("elev")

elev.default <- function(x)
{
    call <- match.call()
    FUN <- function(limit)
        colMeans(sapply(limit, pmin, x = x))
    assign("x", x, env = environment(FUN))
    assign("call", call, env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    FUN
}

elev.grouped.data <- function(x)
{
    call <- match.call()
    FUN <- function(limit)
    {
        ## class in which the limit is situated
        cl <- cut(limit, cj, right = TRUE, label = FALSE)

        ## boundaries of classes strictly below the limit
        xl <- cj[cj <= limit]
        ## boudary of the class including the limit
        xw <- cj[
        ml <- match(x, cj)
        mu <- match(setdiff(cj, x), cj)
        nu <- sum(nj[mu])
        n <- sum(nj[l]) + nu
        (drop(crossprod(x[-length(x)] + diff(x)/2, nj[ml])) +
         limit * nu)/n
    }
    assign("cj", x$cj, env = environment(FUN))
    assign("nj", x$nj, env = environment(FUN))
    assign("call", call, env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    FUN
}
