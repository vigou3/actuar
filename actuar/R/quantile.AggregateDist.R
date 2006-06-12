quantile.AggregateDist <- function(x, approx.lin = TRUE,
                                   p = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
                                   names = TRUE, ...)
{
    n <- length(x)
    upper <- sapply(p, function(q) min(which(x$Fs > q)))
    
    ans <- upper
    if (approx.lin)
    {
        lower <- sapply(p, function(q) max(which(x$Fs < q)))
        h <- (x$Fs[upper] - p) / (x$Fs[upper] - x$Fs[lower])
        ans <- (1 - h) * lower + h * upper
    }
    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(ans) <- formatC(paste(100 * p, "%"), format = "fg", wid = 1, digits = dig)
    }
    ans  
}




