plot.aggregateDist <- function(x, xlim, ...)
{
    ## Function plot() is used for the discretized
    ## CDFs and function curve() in the continuous cases.
    main <- "Aggregate Claims Distribution"
    
    if ("stepfun" %in% class(x)){

        ## Method for class 'ecdf' will most
        ## probably be used.
        NextMethod(ylab = "", main = main) 
    }
    else
    {
        ## Limits for the x-axis are built if none are given
        ## in argument.
        
        if (missing(xlim)){
            mean <- get("mean", environment(x))
            sd <- sqrt(get("var", environment(x)))
            xlim <- c(mean - 3*sd, mean + 3*sd)
        }
        curve(x, main = main, ylab = "", xlim = xlim, ylim = c(0,1))
    }
    mtext(expression(F*scriptstyle(s)(x)), side = 2, line = 2)
    mtext(get("label", environment(x)), line = 0)
    #mtext(comment(x), line = 0)    
}
