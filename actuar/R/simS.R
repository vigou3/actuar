simS <- function(n, model.freq, model.sev)
{ 
    call <- match.call()
    
    ## Get the frequency simulation function.
    rfreq <- match.fun(paste("r", model.freq$dist, sep = ""))
    ## Set the parameters of the frequency distribution
    formals(rfreq)[names(model.freq$par)] <- lapply(model.freq$par, eval.parent)
    ## Simulation of the number of claims
    N <- rfreq(n)

    
    ## Get the severity simulation function.
    rsev <- match.fun(paste("r", model.sev$dist, sep=""))
    ## Set the parameters of the severity distribution.
    formals(rsev)[names(model.sev$par)] <- lapply(model.sev$par, eval.parent)
    ## Simulation of claim amounts
    x <- sapply(sapply(N, rsev),sum)

    
    FUN <- ecdf(x)  ## Computing the empirical CDF
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("call", call, env = environment(FUN))
    #assign("label", "Approximation by simulation", environment(FUN))
    comment(FUN) <- "Approximation by simulation"
    FUN
}


    
