### ===== actuar: an R package for Actuarial Science =====
###
### Display all values of a matrix of vectors by 'unrolling' the
### object vertically or horizontally.
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

### New generic
severity <- function(x, ...) UseMethod("severity")

### Default method. Currently identical to 'unroll' by lack of a
### better alternative.
severity.default <- function(x, bycol = FALSE, drop = TRUE, ...)
    unroll(x, bycol, drop)
