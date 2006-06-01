frequency.simpf <- function(x, ...)
{
    
    contracts <- seq(length = dim(x$data)[1])
    years <- seq(length = dim(x$data)[2])
    names <- list(paste("   Contract", format(contracts)),
                  paste("Year", format(years)))
    print(matrix(sapply(x$data, length),
                 length(contracts),
                 length(years),
                 dimnames = names))
    invisible(x)
}

