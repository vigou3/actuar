### ===== actuar: an R package for Actuarial Science =====
###
### 'mean' method for 'aggregateDist' objects 
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot
mean.aggregateDist <- function(x, ...)
{
    label <- comment(x)	

    ## Simply return the value of the true mean
    ## given in argument in the case of the Normal
    ## and Normal Power approximations.
    
    if (label %in% c("Normal approximation",
                      "Normal Power approximation"))
        return(eval(expression(mean), environment(x)))
    
    else
        return(crossprod(get("x", environment(x)),
                         c(0, diff(eval(expression(y), environment(x)))))[1])
}
