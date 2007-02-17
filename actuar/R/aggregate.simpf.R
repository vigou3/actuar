### ===== actuar: an R package for Actuarial Science =====
###
### Summary statistics of a portfolio
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

aggregate.simpf <- function(x, by = names(x$nodes), FUN = sum, ...)
{
    level.names <- names(x$nodes)       # level names

    ## The most common case should be to aggregate claim amounts by
    ## node. This case is very simple, so we treat it separately.
    if (identical(by, level.names))
        return(cbind(x$classification,
                     array(sapply(x$data, FUN, ...), dim(x$data),
                           dimnames = dimnames(x$data))))

    ## The other possibilities require to split the data in groups as
    ## specified in argument 'by'. If the last level (years) is in
    ## 'by', then the matrix structure must be retained to make the
    ## summaries. Otherwise, it can just be dropped since summaries
    ## will span the years of observation.
    years <- tail(level.names, 1)       # name of last level (years)

    ## Summaries only by last level (years) do not make sense.
    if (identical(by, years))
        stop(paste("'by' must contain at least one level other than '",
                   years, "'", sep = ""))

    ## Convert the sequence of subscripts into factors by pasting the
    ## digits together.
    rows <- setdiff(by, years)          # groups other than years
    s <- x$classification[, rows, drop = FALSE] # subscripts
    f <- apply(s, 1, paste, collapse = "")      # factors
    xx <- split(x$data, f)                      # split data

    ## Function FUN will have to work on lists from now on.
    fun <- function(x, ...) FUN(unlist(x), ...)

    ## Make summaries
    if (years %in% by)
    {
        xx <- lapply(xx, matrix, ncol = ncol(x$data))
        res <- t(sapply(xx, function(x) apply(x, 2, fun, ...)))
        cols <- colnames(x$data)
    }
    else
    {
        res <- sapply(xx, fun, ...)
        cols <- ""
    }

    ## Return results as a matrix
    structure(cbind(s[match(unique(f), f),], res),
              dimnames = list(NULL, c(rows, cols)))
}

frequency.simpf <- function(x, by = names(x$nodes))
{
    FUN <- function(x) if (identical(x, NA)) NA else length(x)
    aggregate(x, by, FUN)
}
