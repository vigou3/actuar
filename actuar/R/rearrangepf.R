rearrangepf <- function(pf)
{
    ## Number of years of observations.
    years <- ncol(pf)

    ## Matrix of the aggregate claim amounts.
    aggregate <- array(dim=dim(pf), sapply(pf, sum))

    ## Matrix of the claim numbers.
    frequencies <- array(dim=dim(pf), sapply(pf, length))

    ## Matrix of the individual claim amounts for the first n - 1
    ## years. Forming this matrix is complicated by the fact that the
    ## number of claims is potentially different for each contract.
    ##
    ## Total number of claims per contract; use 'drop=FALSE' in case
    ## there is only one contract.
    nclaims <- rowSums(frequencies[,1:(years-1), drop=FALSE])

    ## Initialization of the matrix.
    claims <- matrix(NA, nrow(pf), max(nclaims))

    ## Filling of the matrix, contract per contract, only is positions
    ## where there is a claim.
    for (i in 1:nrow(pf))
    {
        if (0 < (nclaimsi <- nclaims[i]))
            claims[i, 1:nclaimsi] <- unlist(pf[i,])[1:nclaimsi]
    }

    ## Matrix of the individual claim amounts for the last
    ## year. Identical to above.
    nclaims <- frequencies[,years]
    claims.last <- matrix(NA, nrow(pf), max(nclaims))
    for (i in 1:nrow(pf))
    {
        if (0 < (nclaimsi <- nclaims[i]))
            claims.last[i, 1:nclaimsi] <- tail(unlist(pf[i,]), nclaimsi)
    }

    list(aggregate=aggregate,
         frequencies=frequencies,
         severities=list(claims=claims, claims.last=claims.last))
}
