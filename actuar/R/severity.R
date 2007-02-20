### ===== actuar: an R package for Actuarial Science =====
###
### Display all values of a matrix of vectors by 'unrolling' the
### object vertically or horizontally.  The method for class 'simpf'
### specifically handles the matrices of data returned by function
### simpf(), where each element is a vector of claim amounts in a
### particular node.
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

### New generic
severity <- function(x, ...) UseMethod("severity")

### Default method for vectors and matrices of vectors (that is, lists
### withour 'dim' attribute or with a 'dim' attribute of length
### 2)
severity.default <- function(x, bycol = FALSE, ...)
{
    if (length(dim(x)) > 2)
        stop("'x' must be a vector or a matrix")

    if (is.null(dim(x)))
        x <- rbind(x)

    fun <- function(x) if (identical(x, NA)) NA else length(x)
    frequencies <- array(sapply(x, fun), dim = dim(x))

    if (bycol)
    {
        lengths <- colSums(frequencies, na.rm = TRUE)
        mat <- matrix(NA, max(lengths), ncol(x), dimnames = dimnames(x))
        for (i in seq_len(ncol(x)))
            if (0 < (lengthi <- lengths[i]))
                mat[seq_len(lengthi), i] <- unlist(x[!is.na(x[, i]), i])
    }
    else
    {
        lengths <- rowSums(frequencies, na.rm = TRUE)
        mat <- matrix(NA, nrow(x), max(lengths),
                      dimnames = list(rownames(x), NULL))
        for (i in seq_len(nrow(x)))
            if (0 < (lengthi <- lengths[i]))
                mat[i, seq_len(lengthi)] <- unlist(x[i, !is.na(x[i, ])])
    }
    drop(mat)
}

## Method for 'simpf' objects. Two differences with the default
## method: 1) data can be grouped according to the levels of the
## portfolio; 2) data can be split by column (e.g. data from all years
## but the last and data from the last year).
severity.simpf <- function(x, by = head(names(x$node), -1),
                           splitcol = NULL, ...)
{
    level.names <- names(x$nodes)       # level names
    ci <- seq_len(ncol(x$data))         # column indexes

    ## Sanity checks
    if (identical(by, level.names))
    {
        warning("nothing to do")
        return(x)
    }
    if (any(!by %in% level.names))
        stop("invalid 'by' specification")

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

    ## Utility function
    fun <- function(x) unlist(x[!is.na(x)])

    ## Split rows according to the 'by' argument.
    s <- x$classification[, by, drop = FALSE]   # subscripts
    f <- apply(s, 1, paste, collapse = "")      # factors
    s <- s[match(unique(f), f),]                # unique subscripts

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
            if (0 < (nc <- ncol(res.last)))/
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
