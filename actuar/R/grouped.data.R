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
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", x, environment(res))
    assign("nj", y, environment(res))
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
            warning("rows will be extracted in sorted order")
            i <- sort(i)
        }
        res <- as.data.frame(NextMethod(x, i))
        class(res) <- c("grouped.data", class(res))
        cj <- get("cj", environment(x))
        nj <- get("nj", environment(x))
        environment(res) <- new.env()
        assign("cj", cj[sort(unique(c(i, i + 1)))], environment(res))
        assign("nj", nj[i], environment(res))
        return(res)
    }

    ## Extraction of classes column case
    if (identical(j, 1))
    {
        cj <- get("cj", environment(x))
        if (missing(i))
            return(cj)
        return(cj[sort(unique(c(i, i + 1)))])
    }

    ## Let other cases (extraction of frequencies column or invalid
    ## arguments) be handled by "[.data.frame".
    NextMethod("[")
}
