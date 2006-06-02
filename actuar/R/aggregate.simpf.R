aggregate.simpf <- function(x, by.contract = TRUE, by.year = TRUE, FUN = sum, ...)
{
    FUN <- match.fun(FUN)
    clabels <- paste("   Contract", format(seq(length = dim(x$data)[1])))
    ylabels <- paste("   Year", format(seq(length = dim(x$data)[2])))
    
    if (by.contract)
    {
        if (by.year)
        {
            y <- matrix(sapply(x$data, FUN),dim(x$data)[1],
                   dim(x$data)[2])
            y[is.nan(y)] <- 0
            dimnames(y) <- list(clabels, ylabels)
        }
        else
        {
            y <- cbind(sapply(sapply(apply(x$data, 1, c),unlist), FUN))
            dimnames(y) <- list(clabels, NULL)
             
        }
    }
    else
        if (by.year)
        {
            y <- rbind(sapply(sapply(apply(x$data, 2, c),unlist), FUN))
            dimnames(y) <- list(NULL, ylabels)
             
        }
        else
        {
            
            y <- FUN(unlist(apply(x$data, 1, c)))
        }
    y
    
    
}



        
    
    
    
