### ===== actuar: an R package for Actuarial Science =====
###
### Extraction of weights method for 'simpf' objects.
###
### AUTHOR:  Vincent Goulet <vincent.goulet@act.ulaval.ca>

weights.simpf <- function(x, ...)
    cbind(x$classification, x$weights)
