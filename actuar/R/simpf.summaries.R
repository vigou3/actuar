### ===== actuar: an R package for Actuarial Science =====
###
### Computing summary statistics and accessing components of a
### portfolio.
###
### AUTHORS: Louis-Philippe Pouliot, Tommy Ouellet, 
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

aggregate.simpf <- function(x, by = names(x$nodes), FUN = sum, ...)
{
    level.names <- names(x$nodes)       # level names
    nlevels <- length(level.names)      # number of levels
    years <- level.names[nlevels]       # name of last level

    ## Match level names in 'by' to those in the model
    by <- match.arg(by, level.names, several.ok = TRUE)

    ## Version of FUN able to work on lists
    fun <- function(x, ...) FUN(unlist(x), ...)

    ## The most common case should be to aggregate claim amounts by
    ## node. This case being very simple, it is treated separately.
    if (identical(by, level.names))
        return(cbind(x$classification,
                     array(sapply(x$data, FUN, ...), dim(x$data),
                           dimnames = dimnames(x$data))))

    ## Summaries only by last level (years) are also simple to handle.
    if (identical(by, years))
        return(apply(x$data, 2, fun, ...))

    ## The other possibilities require to split the data in groups as
    ## specified in argument 'by'. If the last level (years) is in
    ## 'by', then the matrix structure must be retained to make the
    ## summaries. Otherwise, it can just be dropped since summaries
    ## will span the years of observation.
    ##
    ## Convert the sequence of subscripts into factors by pasting the
    ## digits together. It is important *not* to sort the levels in
    ## case the levels in 'by' are not in the same order as in
    ## 'level.names'.
    rows <- setdiff(by, years)          # groups other than years
    s <- x$classification[, rows, drop = FALSE] # subscripts
    f <- apply(s, 1, paste, collapse = "")      # grouping IDs
    f <- factor(f, levels = unique(f))          # factors
    s <- s[match(levels(f), f), , drop = FALSE] # unique subscripts
    xx <- split(x$data, f)                      # split data

    ## Make summaries
    if (years %in% by)
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
              dimnames = list(NULL, c(rows, cols)))
}

frequency.simpf <- function(x, by = names(x$nodes), ...)
{
    freq <- function(x) if (identical(x, NA)) NA else length(x)
    aggregate(x, by, freq)
}

ratios.simpf <- function(object, ...)
{
    if (is.null(weights(object)))
        aggregate(object)
    else
    {
        col <- ncol(object$classification)
        cbind(object$classification, aggregate(object)[, -(1:col)] / weights(object)[, -(1:col)])
    }
}

severity.simpf <- function(x, by = head(names(x$node), -1),
                           splitcol = NULL, ...)
{
    level.names <- names(x$nodes)       # level names
    ci <- seq_len(ncol(x$data))         # column indexes

    ## Match level names in 'by' to those in the model
    by <- match.arg(by, level.names, several.ok = TRUE)

    ## Sanity checks
    if (identical(by, level.names))
    {
        warning("nothing to do")
        return(x)
    }

    ## Convert character 'splitcol' to numeric and then from numeric
    ## or NULL to boolean.
    if (is.character(splitcol))
        splitcol <- pmatch(splitcol, colnames(x$data), duplicates.ok = TRUE)
    if (is.numeric(splitcol) || is.null(splitcol))
        splitcol <- ci %in% splitcol

    ## Unroll claim amounts by column; simplest case
    if (tail(level.names, 1) %in% by)
    {
        if (length(by) > 1)
            stop("invalid 'by' specification")
        x <- x$data
        res <- NextMethod(bycol = TRUE, drop = FALSE)
        return(list(first = res[, !splitcol],
                    last = if (all(!splitcol)) NULL else res[, splitcol]))
    }

    ## Unrolling per row (or group of rows) is more work. It requires
    ## to split the columns of the matrix first, and then to apply the
    ## unrolling procedure twice (if 'splitcol' != NULL).
    ##
    ## Utility function
    fun <- function(x) unlist(x[!is.na(x)])

    ## Split rows according to the 'by' argument.
    s <- x$classification[, by, drop = FALSE]   # subscripts
    f <- apply(s, 1, paste, collapse = "")      # grouping IDs
    f <- factor(f, levels = unique(f))          # factors
    s <- s[match(levels(f), f), , drop = FALSE] # unique subscripts

    ## Keep the 'splitcol' columns for later use.
    x.last <- x$data[, splitcol]

    ## Unroll the "main" block of columns.
    if (all(splitcol))
        res.first <- NULL
    else
    {
        x <- cbind(lapply(split(x$data[, !splitcol], f), fun))
        res.first <- NextMethod(bycol = FALSE, drop = FALSE)
        res.first <-
            if (0 < (nc <- ncol(res.first)))
            {
                dimnames(res.first) <-
                    list(NULL, paste("claim", seq_len(nc), sep = "."))
                cbind(s, res.first)
            }
            else
                NULL
    }

    ## Unroll the 'splitcol' block of columns.
    if (all(!splitcol))
        res.last <- NULL
    else
    {
        x <- cbind(lapply(split(x.last, f), fun))     # split data
        res.last <- NextMethod(bycol = FALSE, drop = FALSE)
        res.last <-
            if (0 < (nc <- ncol(res.last)))
            {
                dimnames(res.last) <-
                    list(NULL, paste("claim", seq_len(nc), sep = "."))
                cbind(s, res.last)
            }
            else
                NULL
    }

    ## Return the result as a list.
    list(first = res.first, last = res.last)
}

weights.simpf <- function(object, ...)
{
    if (is.null(object$weights))
        NULL
    else
        cbind(object$classification, object$weights)
}
