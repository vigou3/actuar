normal <- function(x, mean, var)
{
    call <- match.call()
    Fs <- pnorm((x - mean)/sqrt(var))
    res <- list(fs = c(0, diff(Fs)), Fs = Fs, X = 1:length(Fs), call = call, FUN = approxfun(Fs))
    class(res) <- "AggregateDist"
    res
}
