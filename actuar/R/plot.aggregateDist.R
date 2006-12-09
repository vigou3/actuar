### ===== actuar: an R package for Actuarial Science =====
###
### 'plot' method for 'aggregateDist' objects
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot

plot.aggregateDist <- function(x, xlim, ...)
{
    ## Function plot() is used for the step cdfs and function curve()
    ## in the continuous cases.
    main <- "Aggregate Claim Amount Distribution"
    ylab <- expression(F[S](x))

    if ("stepfun" %in% class(x))
    {
        ## Method for class 'ecdf' will most probably be used.
        NextMethod(main = main, ylab = ylab, ...)
    }
    else
    {
        ## Limits for the x-axis are supplied if none are given
        ## in argument.
        if (missing(xlim))
        {
            mean <- get("mean", environment(x))
            sd <- sqrt(get("var", environment(x)))
            xlim <- c(mean - 3 * sd, mean + 3 * sd)
        }
        curve(x, main = main, ylab = ylab, xlim = xlim, ylim = c(0, 1), ...)
    }
    mtext(comment(x), line = 0.4)
}
