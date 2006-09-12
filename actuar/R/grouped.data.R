### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon


grouped.data <- function(x, y = NULL, ...)
{
    ## 'data.frame' must contain boundaries in first column and number of
    ## data by class in second column.
    if (class(x) == "data.frame")
    {
        y <- x[, 2]
        x <- x[, 1]
    }
    nx <- length(x)
    ny <- length(y)

    ## First data in 'y' won't be used.
    if (nx - ny > 1 || nx - ny < 0)
        stop("length(x) incorrect")
    if (nx - ny == 1)
        y = c(0, y)

    ## Create an object of class 'grouped.data'.
    res = list(cj = x, nj = y)
    class(res) <- c("grouped.data", class(res))
    attr(res, "call") <- sys.call()
    attr(res, "j") <- FALSE
    res
}

print.grouped.data <- function(x, ...)
{
    ## To formate numbers.
    numform <- function(x) formatC(x, digits = 2, format = "e")
    numformy <- function(x) formatC(x)

    ## Use object created by 'grouped' function
    if (attr(x, "j"))
    {
        x <- x$cj
        x1 <- length(x)
        cat("          cj     ", "\n", paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "\n", sep = ""))
    }
    else
    {
        y <- x$nj
        x <- x$cj
        x1 <- length(x)
        cat("          cj     ", "          nj       ", "\n",paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "       ", numformy(y[-1]), "\n", sep = ""))
    }
}




