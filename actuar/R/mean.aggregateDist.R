mean.aggregateDist <- function(x, ...)
{
    #label <- get("label", environment(x))
    label <- comment(x)	

    ## Simply return the value of the true mean
    ## given in argument in the case of the Normal
    ## and Normal Power approximations.
    
    if (label %in% c("Normal approximation",
                      "Normal Power approximation"))
        return(get("mean", environment(x)))
    
    else
        return(crossprod(get("x", environment(x)),
                         c(0, diff(get("y", environment(x))))))
}
        

