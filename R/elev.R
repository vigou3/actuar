### ===== actuar: an R package for Actuarial Science =====
###
### Sample empirical limited value functions for individual and
### grouped data.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca> and
###          Mathieu Pigeon

elev <- function(x, ...)
{
    Call <- match.call()
    UseMethod("elev")
}

elev.default <- function(x, ...)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    FUN <- function(limit)
        colMeans(sapply(limit, pmin, x = x))
    environment(FUN) <- new.env()
    assign("x", sort(x), env = environment(FUN))
    assign("n", length(unique(x)), env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- FALSE
    FUN
}

### This function assumes right-closed intervals, but the numerical
### values are identical for left-closed intervals.
elev.grouped.data <- function(x, ...)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    FUN <- function(limit)
    {
        ## Number of classes.
        r <- length(nj)

        ## This is to avoid numerical problems.
        limit <-  pmin(limit, cj[r + 1])

        ## Class in which the limit is located.
        cl <- findInterval(limit, cj, all.inside = TRUE)

        ## Means for all classes below each limit.
        cjt <- head(cj, max(cl))        # upper bounds
        res1 <- sapply(cl - 1, function(n, x)
                       drop(crossprod(head(x, n), head(nj, n))),
                       (head(cjt, -1) + tail(cjt, -1))/2)

        ## Means for classes with each limit.
        cjt <- cj[cl]                   # lower bounds
        njt <- nj[cl]                   # frequencies
        p <- (limit - cjt) / (cj[cl + 1] - cjt) # prop. to take
        res2 <- njt * p * (cjt + limit)/2 + njt * (1 - p) * limit

        ## Means for classes above each limit.
        res3 <- limit * sapply(r - cl, function(n, x) sum(tail(x, n)),
                               tail(nj, -min(cl)))

        ## Total
        (res1 + res2 + res3)/sum(nj)
    }

    environment(FUN) <- new.env()
    assign("cj", eval(expression(cj), env = environment(x)),
           env = environment(FUN))
    assign("nj", x[, 2], env = environment(FUN))
    assign("n", nrow(x), env = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- TRUE
    FUN
}

### Essentially identical to stats::print.ecdf().
print.elev <- function(x, digits = getOption("digits") - 2, ...)
{
    ## Utility function
    numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")

    ## The rest is adapted from ecdf()
    varname <- if (attr(x, "grouped")) "cj" else "x"
    cat("Empirical LEV \nCall: ")
    print(attr(x, "call"), ...)
    n <- length(xx <- eval(parse(text = varname), env = environment(x)))
    i1 <- 1:min(3, n)
    i2 <- if (n >= 4) max(4, n - 1):n else integer(0)
    cat(" ", varname, "[1:", n, "] = ", numform(xx[i1]), if (n > 3) ", ",
        if (n > 5) " ..., ",  numform(xx[i2]), "\n", sep = "")
    invisible(x)
}

### Essentially identical to stats::summary.ecdf().
summary.elev <- function (object, ...)
{
    cat("Empirical LEV:\t ", eval(expression(n), env = environment(object)),
        "unique values with summary\n")
    summary(knots(object), ...)
}

### Identical to stats::knots.stepfun().
knots.elev <- function(Fn, ...)
{
    if (attr(Fn, "grouped"))
        eval(expression(cj), env = environment(Fn))
    else
        eval(expression(x), env = environment(Fn))
}

plot.elev <- function(x, ..., main = NULL, xlab = "x", ylab = "Empirical LEV")
{
    if (missing(main))
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) cl else sys.call())
        }

    kn <- knots(x)
    Fn <- x(kn)
    plot(kn, Fn,  ..., main = main, xlab = xlab, ylab = ylab)
}
