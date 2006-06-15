simS <- function(x, model.freq, model.sev)
{
    call <- match.call()
    rfreq <- match.fun(paste("r", model.freq$dist, sep = ""))
    formals(rfreq)[names(model.freq$par)] <- lapply(model.freq$par, eval.parent)
    N <- rfreq(x)

    rsev <- match.fun(paste("r", model.sev$dist, sep=""))
    formals(rsev)[names(model.sev$par)] <- model.sev$par
    X <- sapply(sapply(N, rsev),sum)
    X <- sort(X)
    Fs <- ecdf(X)(seq(0, quantile(X, 0.999), by = h))
    fs <- c(0, diff(Fx))
    res <- list(fs = fs, Fs = Fs, call = call, FUN = approxfun(Fs))
    class(res) <- "AggregateDist"
    res
}


    
