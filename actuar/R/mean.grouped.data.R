### ===== actuar: an R package for Actuarial Science =====
###
### Mean (TODO: and summaries) of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

mean.grouped.data <- function(x, ...)
{
    ## Sanity check
    if (!inherits(x, "grouped.data"))
        stop("wrong method")

    cj <- eval(expression(cj), env = environment(x))

    midpoints <- cj[-length(cj)] + diff(cj)/2

    sapply(x[-1], function(x) drop(crossprod(x, midpoints))/sum(x))
}
