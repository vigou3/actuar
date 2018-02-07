### actuar: Actuarial Functions and Heavy Tailed Distributions
###
### Creation of grouped data objects.
###
### The function can create a grouped data object from two types of
### arguments.
###
### 1. Individual data. The call has at least two elements in '...'.
###    The first is then the vector of group boundaries and the others
###    are vectors (or a matrix) of group frequencies.
###
### 2. Group boundaries and frequencies. The call has one or more
###    elements in '...' and either 'breaks' or 'nclass' is provided
###    or 'group' is TRUE. In this case, elements of '...' are grouped
###    using graphics:::hist automatically based on the first element
###    of '...', or with group boundaries 'breaks' if the latter is a
###    vector.
###
### For details on grouped data, see Klugman, Panjer & Willmot, Loss
### Models, Wiley, 1998.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon, Louis-Philippe Pouliot
###
### CREDITS: Manipulation and creation of names taken in part from R
### function data.frame(). Arguments, 'breaks', 'nclass' and their
### treatment taken from R function hist().

grouped.data <- function(..., breaks = "Sturges",
                         include.lowest = TRUE, right = TRUE,
                         nclass = NULL, group = FALSE,
                         row.names = NULL, check.rows = FALSE,
                         check.names = TRUE)
{
    ## Utility function to format numbers.
    numform <- function(x, w)
        formatC(x, digits = 2, width = w, format = "fg")

    ## Keep the calls in '...' in object 'ox', and the evaluated
    ## elements in '...' in object 'x'.
    ox <- as.list(substitute(list(...)))[-1L]
    x <- list(...)
    xlen <- length(x)                   # number of arguments in '...'
    use.br <- !missing(breaks)          # 'breaks' specified

    ## If any elements of '...' are unnamed, set names based on the
    ## variable name provided in the function call (e.g. f(x) -> "x")
    ## or from the deparsed expression (e.g. f(1:3) -> "1:3").
    xnames <- names(x)
    if(length(xnames) != xlen)
	xnames <- character(xlen)
    no.xn <- !nzchar(xnames)
    if (any(no.xn))
    {
        for (i in which(no.xn))
            xnames[i] <- deparse(ox[[i]], nlines = 1L)[1L]
        names(x) <- xnames
    }

    ## Single argument implies individual data.
    if (xlen == 1L)
        group <- TRUE

    ## Avoid using calling 'hist' with 'nclass' specified.
    if (use.br)
    {
        if (!missing(nclass))
            warning("'nclass' not used when 'breaks' is specified")
        if (!(missing(group) || group))
            warning("'group' ignored when 'breaks' is specified")
        group <- TRUE
    }
    else if (!is.null(nclass) && length(nclass) == 1L)
    {
        breaks <- nclass
        if (!(missing(group) || group))
            warning("'group' ignored when 'nclass' is specified")
        group <- TRUE
    }

    if (group)        # individual data in argument; group with 'hist'
    {
        ## Set group boudaries (and the first set of group
        ## frequencies) using the first argument in '...'.
        y <- hist(x[[1]], plot = FALSE, breaks = breaks,
                  include.lowest = include.lowest,
                  right = right)
        br <- y$breaks
        y <- y$counts

        ## If there are other vectors in '...', compute group
        ## frequencies using 'hist' with the group boundaries
        ## determined above. If 'breaks' were set automatically, there
        ## is a great risk of error, but we can't do much better.
        if (xlen > 1)
        {
            f <- function(x, br)
                hist(x, plot = FALSE, breaks = br,
                     include.lowest = include.lowest,
                     right = right)$counts
            y <- cbind(y, sapply(x[-1], f, br = br))
        }

        y <- as.data.frame(y)
        x <- as.data.frame(br)
        names(y) <- xnames
        xnames <- ""
        nx <- nrow(x)
    }
    else              # group boundaries and frequencies in argument
    {
        y <- as.data.frame(x[-1L])      # group frequencies
        x <- as.data.frame(x[[1L]])     # group boundaries
        nx <- nrow(x)

        ## There must be exactly one more group boundary than frequencies.
        if (nx - nrow(y) != 1L)
            stop("invalid number of group boundaries and frequencies")

        ## Replace missing frequencies by zeros.
        nax <- is.na(x)
        if (any(nax))
        {
            x[nax] <- 0
            warning("missing frequencies replaced by zeros")
        }
    }

    ## Return a data frame with formatted group boundaries in the
    ## first column.
    w <- max(nchar(x[-1L, ]))            # longest upper boundary
    xfmt <- paste(if (right) "(" else "[",
                  numform(x[-nx, ], -1), ", ", numform(x[-1L, ], w),
                  if (right) "]" else ")",
                  sep = "")
    res <- data.frame(xfmt, y, row.names = row.names, check.rows = check.rows,
                      check.names = check.names)
    names(res) <- c(xnames[1L], names(y))
    class(res) <- c("grouped.data", "data.frame")
    environment(res) <- new.env()
    assign("cj", unlist(x, use.names = FALSE), environment(res))
    attr(res, "right") <- right
    res
}
