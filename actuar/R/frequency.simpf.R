frequency.simpf <- function(x, by = c("contract", "year"), ...)
{
    MARGIN <- match(by, c("contract", "year"))
    apply(x$data, MARGIN, function(x) length(unlist(x)))
}


