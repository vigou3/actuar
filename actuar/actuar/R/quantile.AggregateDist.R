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

plot.AggregateDist <- function(x)
{
    
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    plot.default(x$fs, main = "Discretized Density Function", xlab = "S", ylab = "f_s(x)", type = "h")
    plot.default(x$Fs, main = "Empirical Cumulative Distribution Function",
                 xlab = "S", ylab = "F_s(x)", type = "s")
    #plot.default(x$fs, main = "Discretized Density Function", xlab = "S", ylab = "f_s(x)", type = "l")
    #plot.default(x$Fs, main = "Empirical Cumulative Distribution Function",
                 #xlab = "S", ylab = "F_s(x)", type = "l", add = TRUE)
    
    plot.default(x$Fs, main = "Empirical Cumulative Distribution Function",
                 xlab = "S", ylab = "F_s(x)", type = "l", xlim = c(0,quantile(x, p=0.9999)))
    lines(diff(c(0,x$Fs)), col = "blue")
    abline(h=0, col = "grey")
    layout(1)
    
}





