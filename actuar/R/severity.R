### ===== actuar: an R package for Actuarial Science =====
###
### Display all values of a matrix of vectors by 'unrolling' the
### object vertically or horizontally.
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

### New generic
severity <- function(x, ...) UseMethod("severity")

### Default method for vectors and matrices of vectors (that is, lists
### without 'dim' attribute or with a 'dim' attribute of length
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
    mat
}
