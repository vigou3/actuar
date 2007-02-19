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

### Create a new generic
severity <- function(x, ...) UseMethod("severity")

### Default method new matrix is created where all the data is
### displayed and empty spaces are filled with NA.
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
        mat <- matrix(NA, max(lengths), ncol(x))
        for (i in seq_len(ncol(x)))
            if (0 < (lengthi <- lengths[i]))
                mat[seq_len(lengthi), i] <- unlist(x[!is.na(x[, i]), i])
    }
    else
    {
        lengths <- rowSums(frequencies, na.rm = TRUE)
        mat <- matrix(NA, nrow(x), max(lengths))
        for (i in seq_len(nrow(x)))
            if (0 < (lengthi <- lengths[i]))
                mat[i, seq_len(lengthi)] <- unlist(x[i, !is.na(x[i, ])])
    }
    drop(mat)
}

## 'simpf' objects are treated more specifically, being displayed more
## appropriately with respect to the actuarial context. The
## possibility of treating separetly the last years of experience is
## also given.

severity.simpf <- function(x, by = head(x$node, -1), split = NULL, ...)
{
    if (is.numeric(split))
        if (any(split < 0 | split > ncol(x$data)))
            stop("invalid split column specification")
    if (is.character(split))
        if (any(is.na(split <- match(split, colnames(x$data)))))
            stop("invalid split column specification")

    level.names <- names(x$nodes)       # level names

    if (identical(by, level.names))
    {
        warning("nothing to do")
        return(x)
    }

    if (tail(level.names, 1) %in% by)
    {
        if (length(by) > 1)
            stop(paste("'", tail(level.names, 1), "' cannot be combined with another level"))
        res <- NextMethod(object = x$data, bycol = TRUE)
    }
    else
    {
        s <- x$classification[, by, drop = FALSE]   # subscripts
        f <- apply(s, 1, paste, collapse = "")      # factors
        xx <- cbind(lapply(split(x$data, f), unlist)) # split data
        res <- NextMethod(object = xx)
    }
    res
}

    x <- x$data
    y <- x
    if (y.exclude %% 1 != 0 | y.exclude < 0)
        stop("y.exclude must be a positive integer")
    if (!bycol) names.first <- list(unlist(dimnames(x)[1]), NULL)
    else names.first <- list(NULL, unlist(dimnames(x)[2])[1:(ncol(x) - y.exclude)])
    x <- x[ , seq(length = ncol(x) - y.exclude)]
    res <- NextMethod()
    dimnames(res) <- names.first

    if (y.exclude > 0)
    {
        if (!bycol) names.last <- names.first
        else names.last <- list(NULL, unlist(dimnames(y)[2])[(ncol(y) - y.exclude + 1):(ncol(y))])

        x <- y[ , seq(ncol(y) - y.exclude + 1, ncol(y))]
        if (!is.matrix(x)) x <- as.matrix(x)
        mat.last <- NextMethod()
        dimnames(mat.last) <- names.last
        res <- list(mat.first = res, mat.last = mat.last)
    }
    print(res)

}
