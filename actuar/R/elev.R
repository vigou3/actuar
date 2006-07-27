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
        ## class in which the limit is located
        cl <- cut(limit, cj, right = TRUE, labels = FALSE)

        ## means for all classes below the limit
        cjt <- head(cj, cl)             # up to the upper limit
        # valid in R 2.4.0 only
        # res1 <- drop(crossprod(((head(cjt, -1) + tail(cjt, -1))/2,
        #                         head(nj, cl - 1))
        res1 <- drop(crossprod((cjt[-1] + cjt[-length(cjt)])/2,
                               head(nj, cl - 1)))

        ## mean for class with the limit
        cjt <- cj[cl]                   # lower limit of the class
        njt <- nj[cl]                   # frequency in the class
        p <- (limit - cjt) / (cj[cl + 1] - cjt) # prop. on each side
        res2 <- njt * p * mean(c(cjt, limit)) + njt * (1 - p) * limit

        ## means for all classes above the limit
        # valid with r 2.4.0 only
        # res3 <- sum(limit * tail(nj, -cl))
        res3 <- sum(limit * nj[-seq(to = cl)])

        ## Total
        (res1 + res2 + res3)/sum(nj)
    }
    assign("cj", x$cj, env = environment(FUN))
    assign("nj", x$nj[-1], env = environment(FUN))
    assign("call", call, env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    FUN
}
