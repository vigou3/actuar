### ===== actuar: an R package for Actuarial Science =====
###
### Creation and manipulation of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon,
### Louis-Philippe Pouliot

grouped.data <- function(..., row.names = NULL, check.rows = FALSE,
                check.names = TRUE)
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
    xnames <- names(x)

    y <- as.data.frame(x[-1])
    x <- as.data.frame(x[[1]])

    nx <- nrow(x)
    ny <- nrow(y)

    ## There must be exactly one class boudary more than frequencies.
    if (nx - ny != 1)
        stop("incorrect number of class boundaries and frequencies")

    ## Return a data frame with formatted class boundaries in the
    ## first column.
    xfmt <- paste("[", numform(x[-nx]), ", ", numform(x[-1]), ")", sep = "")
    res <- data.frame(xfmt, y, row.names = row.names, check.rows = check.rows,
                      check.names = check.names)
    names(res) <- c(xnames[1], names(y))
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", x, environment(res))
    res
}

"[.grouped.data" <- function(x, i, j)
{
    ## Only columns to extract are specified.
    if (nargs() < 3)
    {
        if (missing(i))
            return(x)
        if (is.matrix(i))
            return(as.matrix(x)[i])
        if (identical(1, seq(ncol(x))[i]))
            return(eval(expression(cj), env = environment(x)))
        return(NextMethod("["))
    }

    ## We need row and column indexes to be strictly positive integers.
    ii <- seq(nrow(x))                  # default: all rows
    ii <- if (!missing(i)) ii[i]        # rows specified
    ij <- if (missing(j)) integer(0) else seq(ncol(x))[j] # columns

    ## Extraction of at least the class boundaries (the complicated case).
    if (1 %in% ij)
    {
        ## Extraction of class boundaries in increasing order only
        ## (untractable otherwise).
        if (is.unsorted(ii))
        {
            warning("rows extracted in increasing order")
            ii <- sort(ii)
        }

        ## Fetch the appropriate class boundaries.
        cj <- eval(expression(cj), env = environment(x))
        cj <- cj[sort(unique(c(ii, ii + 1)))]

        ## Extraction of the first column only: return the vector of class
        ## boundaries.
        if (identical(ij, 1))
            return(cj)

        ## Return a modified 'grouped.data' object.
        res <- as.data.frame(NextMethod("["))
        class(res) <- c("grouped.data", class(res))
        environment(res) <- new.env()
        assign("cj", cj, environment(res))
        return(res)
    }

    ## All other cases handled like a regular data frame.
    NextMethod("[")
}

"[<-.grouped.data" <- function(x, i, j, value)
{
    ## We need row and column indexes to be strictly positive integers.
    ii <- seq(nrow(x))
    if (!missing(i))
    {
        ## Matrix indexing only supported for logical matrices.
        if (is.logical(i) && is.matrix(i) && all(dim(i) == dim(x)))
        {
            j <- apply(i, 2, any)       # columns with replacements
            i <- i[, j]                 # keep appropriate column only
        }
        if (is.matrix(i))
            stop("only logical matrix subscripts are allowed in replacement")
        ii <- ii[i]                     # conversion
    }
    ij <- if (missing(j)) integer(0) else c(1, 2)[j]

    ## Replacement in both columns at the same time not supported.
    if (!length(ij) || identical(sort(ij), c(1, 2)))
        stop("impossible to replace boundaries and frequencies simultaneously")

    ## Replacement of class boundaries
    if (identical(ij, 1))
    {
        cj <- eval(expression(cj), env = environment(x))
        cj[sort(unique(c(ii, ii + 1)))] <- value

                                        #ni <- length(i)
                                        #if (length(x) - ni != 1)
                                        #    stop("incorrect number of class boundaries")
        res <- grouped.data(cj, x[, 2])
        names(res) <- names(x)
        return(res)
    }

    ## All other cases handled as a regular data frame.
    NextMethod("[<-")
}
