
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

grData <- function(x, y = NULL, ...)
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

    numform <- function(x) formatC(x, digits = 2, format = "e")
    numformy <- function(x) formatC(x)

    res <- data.frame(Interval = paste("[",numform(x[-length(x)]),", ",numform(x[-1]),")", sep = ""),
                           Frequency = paste(numformy(y[-1]), sep = ""))
    ## Create an object of class 'grouped.data'.
    #res = list(cj = x, nj = y)
    class(res) <- c("grData", "data.frame")
    attr(res, "call") <- sys.call()
    environment(res) <- new.env()
    assign("x", x, envir = environment(res))
    assign("y", y, envir = environment(res))
    res
}



