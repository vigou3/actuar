### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Standard deviation (TODO: and summaries) of grouped data objects
###
### See Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
###
### AUTHOR: Walter Garcia-Fontes <walter.garcia@upf.edu>

sd.grouped.data <- function(x, ...)
{
    ## Square root of variance
    drop(sqrt(var.grouped.data(x)))
}
