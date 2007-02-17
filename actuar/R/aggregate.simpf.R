### ===== actuar: an R package for Actuarial Science =====
###
### Summary statistics of a portfolio
###
### AUTHORS: Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

aggregate.simpf <- function(x, FUN = sum, ...)
{
    cbind(x$classification,
          array(sapply(x$data, FUN, ...), dim(x$data),
                dimnames = dimnames(x$data)))
}

frequency.simpf <- function(x)
{
    FUN <- function(x) if (identical(x, NA)) NA else length(x)
    aggregate(x, FUN)
}


### JUNKYARD
###
### At some point I (VG) wanted aggregate.simpf() to have a 'by'
### argument like the version LPP wrote for the non hierarchical
### version of simpf(). Although feasible, it was getting cumbersome
### and I decided to give up since having a 'by' argument did not seem
### that useful after all. The code I wrote is kept here as archives.
###
### Work needed: come up with a reliable way to create unique IDs from
### a set of subscripts. The scheme explained in the comments is not
### robust should there be more than 10 nodes in a level...

## aggregate.simpf2 <- function(x, by = names(x$nodes), FUN = sum, ...)
## {
##     level.names <- names(x$nodes)       # level names

##     fun <- function(x, ...) FUN(unlist(x), ...) # handle lists

##     ## The most common case should be to aggregate claim amounts by
##     ## node. This case is very simple, so we treat it separately.
##     if (identical(by, level.names))
##         return(cbind(x$classification, apply(x$data, c(1, 2), fun, ...)))

##     ## The other possibilities require to split the data in groups as
##     ## specified in argument 'by'. If the last level (years) is in
##     ## 'by', then the matrix structure must be retained to make the
##     ## summaries. Otherwise, it can just be dropped since summaries
##     ## will span the years of observation.
##     years <- tail(level.names, 1)       # name of last level (years)

##     ## Summaries only by last level (years) do not make sense.
##     if (identical(by, years))
##         stop(paste("'by' must contain at least one level other than '",
##                    years, "'", sep = ""))

##     ## Create unique ids for each grouping specified in 'by'. One row
##     ## of 'x$classification' provides a unique set of subscripts for
##     ## each entity, but this cannot be used directly. The set of
##     ## subscripts, say (i, j, k), is converted into (100i + 10j + k),
##     ## a unique ID of each grouping.
##     rows <- setdiff(by, years)          # groups other than years
##     id <- x$classification[, rows, drop = FALSE] # subscripts
##     f <- drop(crossprod(t(id), 10^((length(rows) - 1):0))) # IDs

##     xx <- split(x$data, f)              # split data

##     if (years %in% by)
##     {
##         xx <- lapply(xx, matrix, ncol = ncol(x$data))
##         res <- t(sapply(xx, function(x) apply(x, MARGIN = 2, FUN = fun, ...)))
##     }
##     else
##         res <- sapply(xx, fun, ...)
##     cbind(y[match(unique(f), f),], res)
## }
