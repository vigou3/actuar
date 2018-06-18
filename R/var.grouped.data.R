### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Variance (TODO: and summaries) of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Walter Garcia-Fontes <walter.garcia@upf.edu>

var.grouped.data <- function(x, ...)
{
    ## Get group boundaries
    cj <- eval(expression(cj), envir = environment(x))

    ## Compute group midpoints
    midpoints <- cj[-length(cj)] + diff(cj)/2

    ## Compute midpoints minus mean and square it
    midsquare <- (midpoints - mean(x))^2

    ## Drop the boundaries column and convert to matrix for use in
    ## crossprod()
    x <- as.matrix(x[-1])

    ## Compute mean per column
    drop(crossprod(x, midsquare))/(colSums(x)-1)
}

