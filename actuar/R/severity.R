### ===== actuar: an R package for Actuarial Science =====
###
### Display all values of a matrix of vectors by 'unrolling' the
### object vertically or horizontally.  The method for class 'simpf'
### specifically handles the matrices of data returned by function
### simpf(), where each element is a vector of losses during a
### particular year.
###
### AUTHORS:  Louis-Philippe Pouliot, Vincent Goulet <vincent.goulet@act.ulaval.ca>

severity <- function(x, bycol = FALSE, ...) UseMethod("severity")

## A new matrix is created where all the data is displayed and empty
## spaces are filled with NA.

severity.default <- function(x, bycol = FALSE, ...)
{
    if (!bycol)
    {
        ncolumns <- ncol(x)
        frequencies <- array(dim=dim(x), sapply(x, length))
        lengths <- rowSums(frequencies[,1:ncolumns, drop=FALSE])
        mat <- matrix(NA, nrow(x), max(lengths))
        for (i in 1:nrow(x))
        {
            if (0 < (lengthi <- lengths[i]))
                mat[i, 1:lengthi] <- unlist(x[i,])[1:lengthi]
        }
    }
    else
    {
        nrows <- nrow(x)
        frequencies <- array(dim=dim(x), sapply(x, length))
        lengths <- colSums(frequencies[1:nrows, , drop=FALSE])
        mat <- matrix(NA, max(lengths), ncol(x))
        for (i in 1:ncol(x))
        {
            if (0 < (lengthi <- lengths[i]))
                mat[1:lengthi, i] <- unlist(x[ , i])[1:lengthi]
        }
    } 
    mat
}

## 'simpf' objects are treated more specifically, being displayed more
## appropriately with respect to the actuarial context. The
## possibility of treating separetly the last years of experience is
## also given.

severity.simpf <- function(x, bycol = FALSE, y.exclude = 0, ...)
{
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

