plot.aggregateDist <- function(x, h = NULL, ...)
{
    if (is.null(h)) h <- ifelse(!is.null(x$h), x$h, 1)
    else h <- h
    
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    min <- ifelse(any(is.na(x$fs)), 0, 0.001)
    xlim <- quantile.aggregateDist(x, p= max(x$Fs, na.rm = TRUE)*c(min, 0.999))
    X <- seq(1, length(x$fs), length = 500)
    plot.default(X, x$fs[X], xlim = xlim,
                 type = "h", main = "Discretized Density Function", xlab = "S", ylab = "f_s(x)")

    
    plot.default(x$Fs, main = "Empirical Cumulative Distribution Function",
                 xlim = xlim, 
                 xlab = "S", ylab = "Fs(x)", type = "s")
    
    plot.default(x$Fs, main = "Distribution of Total Amount of Claims",
                 xlab = "S", ylab = "Fs(x) / fs(x) (blue)", type = "l", xlim = xlim)
    lines(X, 1/h*x$fs[X], col = "blue")
    abline(h=0, col = "grey")
    layout(1) ##################################################################################mod
    
}
