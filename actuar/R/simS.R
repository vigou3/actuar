simS <- function(x, model.freq, model.sev)
{
    rfreq <- match.fun(paste("r", model.freq$dist, sep = ""))
    formals(rfreq)[names(model.freq$par)] <- lapply(model.freq$par, eval.parent)
    N <- rfreq(x)

    rsev <- match.fun(paste("r", model.sev$dist, sep=""))
    formals(rsev)[names(model.sev$par)] <- model.sev$par
    X <- sapply(sapply(N, rsev), sum)
    X
}




    
