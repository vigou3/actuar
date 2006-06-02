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




    
    
    
