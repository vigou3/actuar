### ===== actuar: an R package for Actuarial Science =====
###
### Sample empirical limited value functions for individual and
### grouped data.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

elev <- function(x, ...)
{
    Call <- match.call()
    UseMethod("elev")
}

elev.default <- function(x)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    FUN <- function(limit)
        colMeans(sapply(limit, pmin, x = x))
    assign("x", sort(x), env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- FALSE
    FUN
}

elev.grouped.data <- function(x)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    FUN <- function(limit)
    {
        ## Number of classes.
        r <- length(nj)

        ## Class in which the limit is located.
        cl <- cut(limit, cj, include.lowest = TRUE, labels = FALSE)

        ## Means for all classes below each limit.
        cjt <- head(cj, max(cl))        # upper bounds
        res1 <- sapply(cl - 1, function(n, x)
                       drop(crossprod(head(x, n), head(nj, n))),
                       (cjt[-1] + cjt[-length(cjt)])/2)

        ## Means for classes with each limit.
        cjt <- cj[cl]                   # lower bounds
        njt <- nj[cl]                   # frequencies
        p <- (limit - cjt) / (cj[cl + 1] - cjt) # prop. to take
        res2 <- njt * p * (cjt + limit)/2 + njt * (1 - p) * limit

        ## Means for classes above each limit.
        res3 <- limit * sapply(r - cl, function(n, x) sum(tail(x, n)),
                               nj[-(1:min(cl))])

        ## Total
        (res1 + res2 + res3)/sum(nj)
    }
    assign("cj", x$cj, env = environment(FUN))
    assign("nj", x$nj[-1], env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- TRUE
    FUN
}

knots.elev <- function(fn, ...)
    eval(expression(cj), env = environment(fn))

print.elev <- function(x, digits = getOption("digits") - 2, ...)
{
    numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")

    varname <- if (attr(x, "grouped")) "cj" else "x"
    cat("Empirical LEV \nCall: ")
    print(attr(x, "call"), ...)
    n <- length(xx <- eval(parse(text = varname), env = environment(x)))
    i1 <- 1:min(3, n)
    i2 <- if (n >= 4) max(4, n - 1):n else integer(0)
    cat(" ", varname, "[1:", n, "] = ", numform(xx[i1]), if (n > 3)
        ", ", if (n > 5)
        " ..., ", numform(xx[i2]), "\n", sep = "")
    invisible(x)
}

plot.elev <- function(x, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, col = 1, ...)
  {
    if (attr(x, "grouped"))  
      {
       stop
      }
    else
      {
        xval <- eval(expression(x), env = environment(x))
        plot(xval, x(xval), main = "Empirical Limited Function", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col, type = "o", pch = 20)
      }
  }

