### ===== actuar: an R package for Actuarial Science =====
###
### Creation and manipulation of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon,
### Louis-Philippe Pouliot

grouped.data <- function(...)
{
    ## Utility function
    numform <- function(x) formatC(x, digits = 2, format = "fg")

    ## The function can be called with either one or two arguments.
    ##
    ## One argument: a data frame containing the class boundaries in
    ## the first column and the frequencies in the second column (the
    ## other columns are ignored).
    ##
    ## Two arguments: the first is the vector of class boundaries and
    ## the second the vector of frequencies.
    x <- list(...)

    ## One argument: a data frame
    if (length(x) == 1 & class(x[[1]]) == "data.frame")
    {
        vnames <- names(x[[1]])
        y <- x[[1]][, 2]
        x <- x[[1]][, 1]
    }
    ## Two arguments: two vectors
    else if (length(x) == 2)
    {
        defvnames <- c("Class", "Frequency")
        vnames <- names(x)
        if (is.null(vnames))
            vnames <- defvnames
        else
        {
            mn <- which(names(x) == "")
            vnames[mn] <- defnames[mn]
        }
        y <- x[[2]]
        x <- x[[1]]
    }
    else
        stop("incorrect number of arguments")

    nx <- length(x)
    ny <- length(y)

    ## There must be exactly one class boudary more than frequencies.
    if (nx - ny != 1)
        stop("incorrect number of class boundaries and frequencies")

    ## Return a data frame with formatted class boundaries in the
    ## first column.
    xfmt <- paste("[", numform(x[-nx]), ", ", numform(x[-1]), ")", sep = "")
    res <- data.frame(xfmt, y)
    names(res) <- vnames
    class(res) <- c("gouped.data", "data.frame")
    attr(res, "call") <- sys.call()
    environment(res) <- new.env()
    assign("cj", x, envir = environment(res))
    assign("nj", y, envir = environment(res))
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
