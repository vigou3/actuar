### ===== actuar: an R package for Actuarial Science =====
###
### Summary statistics of a portfolio
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

aggregate.simpf <- function(x, by = c("contract", "year"), FUN = sum, ...)
{
    FUN <- match.fun(FUN)
    MARGIN <- match(by, c("contract", "year"))
    apply(x$data, MARGIN, function(x) FUN(unlist(x)))
}

frequency.simpf <- function(x, by = c("contract", "year"), ...)
{
    MARGIN <- match(by, c("contract", "year"))
    apply(x$data, MARGIN, function(x) length(unlist(x)))
}
