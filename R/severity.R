### ===== actuar: An R Package for Actuarial Science =====
###
### Display all values of a matrix of vectors by 'unrolling' the
### object vertically or horizontally.
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

### New generic
severity <- function(x, ...) UseMethod("severity")

### Default method. Currently identical to 'unroll' by lack of a
### better alternative. This default method is called (at least) by
### severity.portfolio() with NextMethod() and receives additional,
### unused arguments. As such, it should not check the content of
### '...' with chkDots(). [I'm not quite sure whether this is bad
### programming or not.]
severity.default <- function(x, bycol = FALSE, drop = TRUE, ...)
    unroll(x, bycol, drop)
