frequency.simpf <- function(x, by.contract = TRUE, by.year = TRUE, ...)
{
    contracts <- seq(length = dim(x$data)[1])
    years <- seq(length = dim(x$data)[2])
    mat.freq <- matrix(sapply(x$data, length), nrow(x$data), ncol(x$data))
    if (by.contract)
    {
        if (by.year)
        {
            freq <- mat.freq
            dimnames(freq) <- list(paste("   Contract", format(contracts)),
                                   paste("Year", format(years)))
        }
        else
        {
            freq <- cbind(rowSums(mat.freq))
            dimnames(freq) <- list(paste("   Contract", format(contracts)), NULL)             
        }
    }
    else
        if (by.year)
        {
            freq <- rbind(colSums(mat.freq))
            dimnames(freq) <- list(NULL, paste("Year", format(years)))             
        }
        else freq <- sum(mat.freq)            
    freq
}

