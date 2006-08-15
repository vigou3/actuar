
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

## Calculate an empirical distribution function.
ogive <- function(x, y = NULL)
{
    ## Use object created by 'grouped' function.
    if (class(x) == "grouped.data")
    {
        y <- x$nj
        x <- x$cj
    }
    ## An error message is issued if 'x' is empty.
    if (length(x) < 1)
      stop("'x' must have 1 or more non-missing values")

    ##Create an object of class 'ogive'.
    Fnt <- approxfun(x, cumsum(y) / sum(y), yleft = 0, yright = 1, method = "linear", ties = "ordered")
    class(Fnt) <- c("ogive", class(Fnt))
    attr(Fnt, "call") <- sys.call()
    Fnt
}

print.ogive <- function(x, digits = getOption("digits") - 2, ...)
{
    ## To formate numbers.
    numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")

    ## Create a structure for the presentation of empirical distribution function.
    cat("Empirical CDF for grouped data \nCall: ")
    print(attr(x, "call"), ...)
    nc <- length(xxc <- eval(expression(x), env = environment(x)))
    nn <- length(xxn <- eval(expression(y), env = environment(x)))
    i1 <- 1:min(3, nc)
    i2 <- if (nc >= 4) max(4, nc - 1):nc else integer(0)
    i3 <- 1:min(3, nn)
    i4 <- if (nn >= 4) max(4, nn - 1):nn else integer(0)
    cat(" cj = ", numform(xxc[i1]), if (nc > 3) ", ", if (nc > 5) " ..., ", numform(xxc[i2]), "\n", sep = "")
    cat(" nj = ", numform(xxn[i1]), if (nn > 3) ", ", if (nn > 5) " ..., ", numform(xxn[i2]), "\n", sep = "")
    invisible(x)
}

