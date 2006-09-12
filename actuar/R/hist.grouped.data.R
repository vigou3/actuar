
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon



## Calculate an empirical density function and create the histogram.
hist.grouped.data <- function(x, y = NULL, main = "Histogram", xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "f(x)", plot = TRUE, ...)
{
    ## Use object created by 'grouped' function.
    if (inherits(x, "grouped.data")
    {
        y <- x$nj
        x <- x$cj
    }

    ## Create an object of class 'histogram'.
    fnt <- approxfun(x, c(0, y[-1] / (sum(y) * diff(x))), yleft = 0, yright = 0, f = 1, method = "constant")
    r <- structure(list(cj = x, nj = y[-1], density = fnt(x)), class = "histogram")

    ## If 'plot' is true, histogram is created, else, boudaries, number of data by class and density are returned.
    if (plot)
    {
        plot(x , fnt(x), main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "S", frame = FALSE)
        segments(x, 0, x, fnt(x))
        segments(0, 0, max(x), 0)
        invisible(r)
    }
    else
        r
}

