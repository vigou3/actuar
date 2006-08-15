
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

"[.grouped.data" <- function(x, i, j)
{
    ## If 'i' is missing, all rows are presented.
    if (missing(i))
        i <- (1:(length(x$cj)-1))

    ## If 'j' is missing, all columns are presented.
    if (missing(j))
      {
        nj <- x$nj[-1]
        nj <- c(0, nj[(i)])
        cj <- x$cj
        cj <- cj[c((i), max(i) + 1)]
        j <- 0
        res = list(cj = cj, nj =  nj)
        class(res) <- c("grouped.data", class(res))
        attr(res, "call") <- sys.call()
        attr(res, "j") <- FALSE
      }

    ## Create an object of class 'grouped.data' or 'numeric'.
    if (j == 1)
    {
        cj <- x$cj
        cj <- cj[c((i), max(i) + 1)]
        res = list(cj = cj, j = 1)
        class(res) <- c("grouped.data", class(res))
        attr(res, "call") <- sys.call()
        attr(res, "j") <- TRUE
    }
    if (j == 2)
    {
        nj <- x$nj[-1]
        nj <- nj[(i)]
        res = nj
        class(res) <- c("numeric")
    }
    res
}
