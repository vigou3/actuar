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
    ## the first column and the class frequencies in the second
    ## column. The first value of the second column is dropped. The
    ## other columns, if any, are ignored.
    ##
    ## Two arguments: the first is the vector of class boundaries and
    ## the second the vector of frequencies. The second vector is one
    ## element shorter than the first.
    x <- list(...)
    if (length(x) == 1 & class(x[[1]]) == "data.frame")
    {
        cnames <- names(x[[1]])
        y <- x[[1]][-1, 2]              # drop dummy frequency
        x <- x[[1]][, 1]
    }
    else if (length(x) == 2)
    {
        cnames.def <- c("Class", "Frequency")
        cnames <- names(x)
        if (is.null(cnames))
            cnames <- cnames.def
        else
        {
            mn <- which(names(x) == "")
            cnames[mn] <- defnames[mn]
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
    names(res) <- cnames
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", x, environment(res))
    res
}

"[.grouped.data" <- function(x, i, j)
{
    ## Extraction of both columns case
    if (missing(j) || identical(sort(j), c(1, 2)))
    {
        if (missing(i))
            return(x)
        if (is.unsorted(i))
        {
            warning("rows are extracted in increasing order")
            i <- sort(i)
        }
        res <- as.data.frame(NextMethod("["))
        class(res) <- c("grouped.data", class(res))
        cj <- get("cj", environment(x))
        environment(res) <- new.env()
        assign("cj", cj[sort(unique(c(i, i + 1)))], environment(res))
        return(res)
    }

    ## Extraction of class boudaries column case
    if (j %in% c(1, -2))
    {
        cj <- get("cj", environment(x))
        if (missing(i))
            return(cj)
        return(cj[sort(unique(c(i, i + 1)))])
    }

    ## Leave other cases (extraction of frequencies column or invalid
    ## arguments) to "[.data.frame"().
    NextMethod("[")
}

"[<-.grouped.data" <- function(x, i, j, value)
{
    ## Impossible to assign both columns at the same time
    if (missing(j) || identical(sort(j), c(1, 2)))
        stop("impossible to replace class boundaries and class frequencies simultaneously")

    ## Replacement of class boundaries
    if (j %in% c(1, -2))
    {
        ni <- length(i)
        if (length(x) - ni != 1)
            stop("invalid number of class boundaries")
        cj <- get("cj", environment(x))
        if (missing(i))
            cj <- value
        else
            cj[sort(unique(c(i, i + 1)))] <- value
        res <- grouped.data(cj, x[, 2])
        names(res) <- names(x)
        return(res)
    }

    ## Leave replacement of frequencies to "[<-.data.frame"().
    NextMethod("[<-")
}
