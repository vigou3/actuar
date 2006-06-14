plot.AggregateDist <- function(x, ...)
{
    
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
    plot.default(x$fs, main = "Discretized Density Function",
                 xlab = "S", ylab = "f_s(x)", type = "h", xlim = c(0,quantile.AggregateDist(x, p=0.999)))
    plot.default(x$Fs, main = "Empirical Cumulative Distribution Function",
                 xlab = "S", ylab = "F_s(x)", type = "s")
    
    plot.default(x$Fs, main = "Distribution of S",
                 xlab = "S", ylab = "F_s(x) / f_s(x) (blue)", type = "l",ylim = c(0, 1))
    lines(diff(c(0,x$Fs)), col = "blue")
    abline(h=0, col = "grey")
    layout(1)
    
}
