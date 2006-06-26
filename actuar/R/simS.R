simS <- function(x, model.freq, model.sev, x.scale = 1)
{
    call <- match.call()
    rfreq <- match.fun(paste("r", model.freq$dist, sep = ""))
    formals(rfreq)[names(model.freq$par)] <- lapply(model.freq$par, eval.parent)
    N <- rfreq(x)

    rsev <- match.fun(paste("r", model.sev$dist, sep=""))
    formals(rsev)[names(model.sev$par)] <- lapply(model.sev$par, eval.parent)
    X <- sapply(sapply(N, rsev),sum)
    X <- x.scale*sort(X)
    FUN <- ecdf(X)
    class(FUN) <- c("aggregateDist", class(FUN))
    assign("call", call, env = environment(FUN))
    assign("X", X, env = environment(FUN))
    assign("x.scale", x.scale, env = environment(FUN))
    FUN
}


    
