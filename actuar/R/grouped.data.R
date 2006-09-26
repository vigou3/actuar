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

    ## The function must be called with at least two arguments. The
    ## first is the vector of class boundaries. The others are vectors
    ## of class frequencies. All arguments will be converted to data
    ## frames.
    x <- list(...)
    xnames <- names(x)                  # preserve names
    y <- as.data.frame(x[-1])           # class frequencies
    x <- as.data.frame(x[[1]])          # class boundaries
    nx <- nrow(x)
    ny <- nrow(y)

    ## There must be exactly one class boundary more than frequencies.
    if (nx - ny != 1)
        stop("incorrect number of class boundaries and frequencies")

    ## Return a data frame with formatted class boundaries in the
    ## first column.
    xfmt <- paste("(", numform(x[-nx, ]), ", ", numform(x[-1, ]), "]",
                  sep = "")
    res <- data.frame(xfmt, y, row.names = row.names, check.rows = check.rows,
                      check.names = check.names)
    names(res) <- c(xnames[1], names(y))
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", unlist(x, use.names = FALSE), environment(res))
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
        res <- NextMethod("[")
        if (identical(seq(ncol(x))[i], as.integer(1)))
            environment(res) <- environment(x)
        return(res)
    }

    ## We need row and column indexes to be strictly positive integers.
    ii <- if (missing(i)) seq(nrow(x)) else seq(nrow(x))[i]
    ij <- if (missing(j)) integer(0) else seq(ncol(x))[j]

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
        if (identical(ij, as.integer(1)))
            return(cj)

        ## Return a modified 'grouped.data' object.
        res <- NextMethod("[")
        #class(res) <- c("grouped.data", class(res))
        environment(res) <- new.env()
        assign("cj", cj, environment(res))
        return(res)
    }

    ## All other cases handled like a regular data frame.
    NextMethod("[")
}

"[<-.grouped.data" <- function(x, i, j, value)
{
    print(nargs())
    nA <- nargs()
    if (nA == 4)
    {
        i <- if (missing(i)) NULL else i
        j <- if (missing(j)) NULL else j
    }
    else if (nA == 3)
    {
        ## No arguments inside [ ]: only replacing by NULL is supported.
        if (missing(i) && missing(j))
        {
            if (is.null(value))
                return(x[logical(0)])
            stop("impossible to replace boundaries and frequencies simultaneously")
        }
        ## Indexing by a logical matrix is supported, but only two
        ## types of replacement are allowed: replacing in the
        ## first column only, or replacing in any column but the
        ## first.
        if (is.logical(i) && is.matrix(i) && all(dim(i) == dim(x)))
        {
            j <- apply(i, 2, any)       # columns with replacements
            if (match(TRUE, j) == 1)    # boundaries to replace
            {
                if (length(j) > 1)      # boundaries and frequencies
                    stop("impossible to replace boundaries and frequencies simultaneously")
                i <- i[, j]             # boundaries only
            }
            return(NextMethod("[<-"))   # frequencies only
        }
        ## Indexing by a non logical matrix is not supported.
        if (is.matrix(i))
            stop("only logical matrix subscripts are allowed in replacement")
        ## Indexing by a vector: the argument specifies columns to
        ## replace.
        j <- i
        i <- NULL
    }
    else
        stop("need 0, 1, or 2 subscripts")

    ## We need row and column indexes to be strictly positive integers.
    ii <- if (is.null(i)) seq(nrow(x)) else seq(nrow(x))[i]
    ij <- if (is.null(j)) integer(0) else seq(ncol(x))[j]

    ## Replacement in the column of boundaries and any other column is
    ## not supported (untractable).
    if (!length(ij) || (length(ij) < 1 && 1 %in% ij))
        stop("impossible to replace boundaries and frequencies simultaneously")

    ## Replacement of class boundaries only.
    if (identical(ij, as.integer(1)))
    {
        cj <- eval(expression(cj), env = environment(x))
        cj[sort(unique(c(ii, ii + 1)))] <- value
        res <- grouped.data(cj, x[, -1])
        names(res) <- names(x)
        return(res)
    }

    ## All other cases handled like a regular data frame.
    NextMethod("[<-")
}
