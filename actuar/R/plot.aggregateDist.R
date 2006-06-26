plot.aggregateDist <- function(x, xlim, ...)
{
    main <- "Distribution of the Aggregate Claims"
    ylab <- "Fs(x)"
    if ("stepfun" %in% class(x)) 
        NextMethod(ylab = ylab, main = main)
    else
    {
        if (missing(xlim)){
            mean <- get("mean", environment(x))
            sd <- sqrt(get("var", environment(x)))
            xlim <- c(mean - 3*sd, mean + 3*sd)
        }
        curve(x, main = main, xlim = xlim)
    }
}
