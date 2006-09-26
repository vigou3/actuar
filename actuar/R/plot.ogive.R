### ===== actuar: an R package for Actuarial Science =====
###
### plot() method for ogives
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

plot.ogive <- function(x, ..., main = NULL, xlab = "x", ylab = "F(x)")
{
    ## Sanity check
    if (!inherits(x, "ogive"))
        stop("wrong method")

    if (missing(main))
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) cl else sys.call())
        }

    kn <- knots(x)
    Fn <- x(kn)
    plot(kn, Fn,  ..., type = "o", pch = 16,
         main = main, xlab = xlab, ylab = ylab)
}
