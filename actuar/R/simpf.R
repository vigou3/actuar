simpf <- function(contracts, years, model.freq, model.sev, weights)
{
    freq <- model.freq
    sev <- model.sev
    ## Assign a matrix of weights if none are given in argument.
    if (missing(weights))
        weights <- matrix(1, contracts, years)

    ## Verify that the dimensions of the weights matrix match the
    ## 'contracts' and 'years' arguments.
    if (!isTRUE(all.equal(dim(weights), c(contracts, years))))
        stop(paste("dimensions of matrix 'weights' should be c(", contracts, ", ", years, ")", sep=""))

    ## Total number of observations in the portfolio (used often).
    nobs <- contracts * years

    ## Simulation of the frequencies. If 'model.freq' is NULL, this is
    ## equivalent to having one claim per contract per
    ## year. Otherwise, the number of claims is simulated for each
    ## contract and each year.
    if (is.null(model.freq))
    {
        N <- rep(1, nobs)
    }
    else
    {
        ## Get the frequency simulation function.
        rfreq <- match.fun(paste("r", model.freq$dist1, sep=""))

        ## The presence of the string "Lambda" in model.freq$par1
        ## indicates a compound model. If present, then we have to
        ## simulate the values of the compounding parameter using the
        ## distribution specified in model.freq$dist2. Otherwise, we
        ## simply simulate from the specified distribution.
        if (is.na(pmatch("Lambda", as.character(model.freq$par1))))
        {
            ## If there is no compounding parameter but a compounding
            ## distribution is specified, issue a warning.
            if (exists("model.freq$dist2"))
                warning("A compounding distribution for the frequency of claims is specified, but no compounding parameter")
        }
        else
        {
            ## Get the compounding distribution simulation function and set
            ## its parameters.
            rlambda <- match.fun(paste("r", model.freq$dist2, sep=""))
            formals(rlambda)[names(model.freq$par2)] <- model.freq$par2

            ## Simulation of the compounding parameters (frequency risk levels).
            Lambda <- rlambda(contracts)
        }


        ## Set the parameters of the frequency distribution by
        ## evaluating the expression in 'Lambda' and/or 'weights'
        ## given in model.freq$par1.
        formals(rfreq)[names(model.freq$par1)] <- lapply(model.freq$par1, eval.parent)

        ## Simulation of the number of claims per year and per
        ## contract (vector Lambda, if any, is recycled).
        N <- rfreq(nobs)
    }

    ## Simulation of the claim amounts. If 'model.sev' is NULL, this
    ## is equivalent to simulating frequencies only. Otherwise, claim
    ## amounts are simulated for each claim.
    if (is.null(model.sev))
    {
        X <- N
    }
    else
    {
        ## Get the severity simulation function.
        rsev <- match.fun(paste("r", model.sev$dist1, sep=""))

        ## The presence of the string "Theta" in model.sev$par1
        ## indicates a compound model. If present, then we have to
        ## simulate the values of the compounding parameter using the
        ## distribution specified in model.sev$dist2. Otherwise, we
        ## simply simulate from the specified distribution.
       if (is.na(pmatch("Theta", as.character(model.sev$par1))))
       {
            ## If there is no compounding parameter but a compounding
            ## distribution is specified, issue a warning.
            if (exists("model.sev$dist2"))
                warning("A compounding distribution for the amount of claims is specified, but no compounding parameter")

            ## Set the parameters of the severity distribution.
            formals(rsev)[names(model.sev$par1)] <- model.sev$par1

            ## Since the parameters of the severity distribution do
            ## not change from one contract to another, we can
            ## immediately simulate all claim amounts.
            X <- sapply(N, rsev)
        }
        else
        {
            ## Get the compounding distribution simulation function and set
            ## its parameters.
            rtheta <- match.fun(paste("r", model.sev$dist2, sep=""))
            formals(rtheta)[names(model.sev$par2)] <- model.sev$par2

            ## Simulation of the compounding parameters (severity risk levels).
            Theta <- rtheta(contracts)

            ## Simulation of claim amounts in the case of a compound
            ## model is more complicated since the severity risk
            ## parameter (potentially) changes from one contract to
            ## another. We must therefore be able to distinguish which
            ## Theta to use for each claim amount simulation. For
            ## this, every occurence of 'Theta' in the severity
            ## distributions parameters must be replaced by 'Theta[i]'.
            model.sev$par1 <- lapply(model.sev$par1, function(x) parse(text=sub("Theta", "Theta[i]", deparse(x))))

            ## We then use an auxiliary function (to be used in
            ## lapply()) to run through all contracts and years. It
            ## will choose the correct Theta to use.
            f <- function(j)
            {
                ## index of the contract to simulate
                i <- 1 + (j - 1) %% contracts

                ## set the parameters of the severity distribution;
                ## this is where the '[i]' pasted above is used
                formals(rsev)[names(model.sev$par1)] <- lapply(model.sev$par1, function(x) {force(i); eval(x)})

                ## simulation of claim amounts for this contract
               rsev(N[j])
            }

            ## Simulation of claim amounts for every contract and year.
            X <- sapply(1:nobs, f)
        }
    }

    ## Return individual claim amounts as a two dimension list or
    ## simple matrix, if possible.
    dim(X) <- c(contracts, years)
    dimnames(X) <- list(paste("Contract", 1:contracts, sep = "."), paste("Year", 1:years, sep = "."))
    res <- list(data = X,
                weights = weights,
                freq = freq,
                sev = sev)
    class(res) <- "simpf"
    res
}

print.simpf <- function(x, ...)
{
    freq1 <- paste(x$freq$dist1, "(",
                   paste(as.character(x$freq$par1), collapse = ", "), ")", sep = "")
    if (!is.null(x$freq$dist2))
        freq2 <- paste("/ ",x$freq$dist2, "(",
                       paste(as.character(x$freq$par2), collapse = ", "), ")", sep = "")
    else freq2 <- NULL

    
    sev1 <- paste(x$sev$dist1, "(", paste(as.character(x$sev$par1), collapse = ", "), ")", sep = "")
    if (!is.null(x$sev$dist2))
        sev2 <- paste("/ ",x$sev$dist2, "(", paste(as.character(x$sev$par2), collapse = ", "), ")", sep = "")
    else sev2 <- NULL

    
    cat("\nPortfolio of claim amounts \n\n")
    cat("  Frequency Model: ", freq1, freq2, "\n")
    cat("  Severity Model: ", sev1, sev2, "\n")
    cat("\n  Number of claims per contract and per year: \n\n")
    print(array(sapply(x$data, length),dim(x$data), dimnames = dimnames(x$data)))
    invisible(x)   
}
