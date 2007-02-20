### ===== actuar: an R package for Actuarial Science =====
###
### Summary statistics of a portfolio
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

aggregate.simpf <- function(x, by = names(x$nodes), FUN = sum, ...)
{
    level.names <- names(x$nodes)       # level names
    nlevels <- length(level.names)      # number of levels

    by <- pmatch(by, level.names, duplicates.ok = TRUE)

    ## Version of FUN able to work on lists
    fun <- function(x, ...) FUN(unlist(x), ...)

    ## The most common case should be to aggregate claim amounts by
    ## node. This case being very simple, it is treated separately.
    if (identical(length(by), nlevels))
        return(cbind(x$classification,
                     array(sapply(x$data, FUN, ...), dim(x$data),
                           dimnames = dimnames(x$data))))

    ## Summaries only by last level (years) are also simple to handle.
    if (identical(by, nlevels))
        return(apply(x$data, 2, fun, ...))

    ## The other possibilities require to split the data in groups as
    ## specified in argument 'by'. If the last level (years) is in
    ## 'by', then the matrix structure must be retained to make the
    ## summaries. Otherwise, it can just be dropped since summaries
    ## will span the years of observation.
    years <- level.names[nlevels]       # name of last level (years)

    ## Convert the sequence of subscripts into factors by pasting the
    ## digits together.
    rows <- setdiff(by, nlevels)        # groups other than years
    s <- x$classification[, rows, drop = FALSE] # subscripts
    f <- apply(s, 1, paste, collapse = "")      # factors
    s <- s[match(unique(f), f),]                # unique subscripts
    xx <- split(x$data, f)                      # split data

    ## Make summaries
    if (nlevels %in% by)
    {
        xx <- lapply(xx, matrix, ncol = ncol(x$data))
        res <- t(sapply(xx, function(x, ...) apply(x, 2, fun, ...), ...))
        cols <- colnames(x$data)
    }
    else
    {
        res <- sapply(xx, fun, ...)
        cols <- deparse(substitute(FUN))
    }

    ## Return results as a matrix
    structure(cbind(s, res),
              dimnames = list(NULL, c(level.names[rows], cols)))
}

frequency.simpf <- function(x, by = names(x$nodes))
{
    FUN <- function(x) if (identical(x, NA)) NA else length(x)
    aggregate(x, by, FUN)
}
