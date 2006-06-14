aggregate.simpf <- function(x, by = c("contract", "year"), FUN = sum, ...)
{
    FUN <- match.fun(FUN)
    MARGIN <- match(by, c("contract", "year"))
    apply(x$data, MARGIN, function(x) FUN(unlist(x)))
}





        
    
    
    
