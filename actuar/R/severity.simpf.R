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


