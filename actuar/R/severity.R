severity <- function(x, ...) UseMethod("severity")


severity.default <- function(x, byrow = TRUE, ...)
{
    if (byrow)
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


severity.simpf <- function(x, byrow = TRUE, y.exclude = 0, ...)
{
    x <- x$data
    y <- x
    
    if (byrow) names.first <- list(paste("Contract", seq(length = nrow(y))), NULL)
    else names.first <- list(NULL, paste("Year", seq(length = ncol(x) - y.exclude)))
    
    x <- x[ , seq(length = ncol(x) - y.exclude)]
    mat.first <- NextMethod()
    dimnames(mat.first) <- names.first
    res <- list(mat.first = mat.first)
    
    if (y.exclude > 0)
    {
        if (byrow) names.last <- names.first
        else names.last <- list(NULL, paste("Year", seq(ncol(y) - y.exclude + 1, ncol(y))))
        
        x <- y[ , seq(ncol(y) - y.exclude + 1, ncol(y))]
        if (!is.matrix(x)) x <- as.matrix(x)
        mat.last <- NextMethod()
        dimnames(mat.last) <- names.last
        res$mat.last <- mat.last       
    }
    res
}
