### ===== actuar: An R Package for Actuarial Science =====
###
### Ogive for grouped data.
###
### A default method exists for either a vector of individual data or
### two vectors of group boundaries and group frequencies. It first
### creates a grouped data object using 'grouped.data' and then
### dispatches to the main method to get the ogive.
###
### For the definition of the ogive, see Klugman, Panjer & Willmot,
### Loss Models, Wiley, 1998.
###
### More details on the admissible arguments for the default method
### are to be found in ./grouped.data.R.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon
###
### CREDITS: Arguments, 'breaks', 'nclass' and their treatment taken
### from R function hist().

ogive <- function(x, ...)
{
    Call <- match.call()
    UseMethod("ogive")
}

ogive.default <- function(x, y = NULL,
                          breaks = "Sturges", nclass = NULL, ...)
{
    chkDots(...)           # method does not use '...'

    ## Avoid using calling 'hist' with 'nclass' specified.
    if (!missing(breaks))
    {
        if (!missing(nclass))
            warning("'nclass' not used when 'breaks' is specified")
    }
    else if (!is.null(nclass) && length(nclass) == 1L)
        breaks <- nclass

    ## Create the "grouped.data" object.
    x <- if (is.null(y))   # one argument: individual data
             grouped.data(x, breaks = breaks)
         else              # two arguments: boundaries and frequencies
             grouped.data(x, y)

    ## Get the ogive from the "grouped.data" method.
    ogive.grouped.data(x)
}

ogive.grouped.data <- function(x, ...)
{
    chkDots(...)                      # method does not use '...'

    ## Group frequencies in the second column of the data frame; group
    ## boundaries in the environment of 'x'.
    y <- x[, 2L]
    x <- eval(expression(cj), envir = environment(x))

    ## Create an object of class 'ogive'.
    res <- approxfun(x, cumsum(c(0, y)) / sum(y), yleft = 0, yright = 1,
                     method = "linear", ties = "ordered")
    class(res) <- c("ogive", class(res))
    attr(res, "call") <- sys.call()
    res
}

### Essentially identical to stats::print.ecdf().
print.ogive <- function(x, digits = getOption("digits") - 2, ...)
{
    ## Utility function
    numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")

    ## The rest is adapted from ecdf()
    cat("Ogive for grouped data \nCall: ")
    print(attr(x, "call"), ...)
    nc <- length(xxc <- get("x", envir = environment(x)))
    nn <- length(xxn <- get("y", envir = environment(x)))
    i1 <- 1L:min(3L, nc)
    i2 <- if (nc >= 4L) max(4L, nc - 1L):nc else integer(0)
    i3 <- 1L:min(3L, nn)
    i4 <- if (nn >= 4L) max(4L, nn - 1L):nn else integer(0)
    cat("    x = ", numform(xxc[i1]), if (nc > 3L) ", ",
        if (nc > 5L) " ..., ", numform(xxc[i2]), "\n", sep = "")
    cat(" F(x) = ", numform(xxn[i3]), if (nn > 3L) ", ",
        if (nn > 5L) " ..., ", numform(xxn[i4]), "\n", sep = "")
    invisible(x)
}

### Essentially identical to stats::summary.ecdf().
summary.ogive <- function (object, ...)
{
    n <- length(eval(expression(x), envir = environment(object)))
    header <- paste("Ogive:	 ", n,
                    "unique values with summary\n")
    structure(summary(knots(object), ...),
              header = header, class = "summary.ogive")
}

### Identical to stats::print.summary.ecdf().
print.summary.ogive <- function(x, ...)
{
    cat(attr(x, "header"))
    y <- unclass(x); attr(y, "header") <- NULL
    print(y, ...)
    invisible(x)
}

### Identical to stats::knots.stepfun().
knots.ogive <- stats:::knots.stepfun

plot.ogive <- function(x, main = NULL, xlab = "x", ylab = "F(x)", ...)
{
    if (missing(main))
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) cl else sys.call())
        }

    kn <- knots(x)
    Fn <- x(kn)
    plot(kn, Fn,  ..., type = "o", pch = 16,
         main = main, xlab = xlab, ylab = ylab)
}
