### ===== actuar: an R package for Actuarial Science =====
###
### Extract individual claim amounts from an object of class
### 'simpf'. Two differences with the default method: 1) data can be
### grouped according to the levels of the portfolio; 2) data can be
### split by column (e.g. data from all years but the last and data
### from the last year).
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

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
        res <- NextMethod(bycol = TRUE)
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
        res.first <- NextMethod(bycol = FALSE)
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
        res.last <- NextMethod(bycol = FALSE)
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
