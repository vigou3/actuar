### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the loss distributions facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar);
if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz"))
            mfrow = c(2, 2))


###
### PROBABILITY LAWS
###

## The package provides "d", "p", "q" and "r" functions for all the
## probability laws useful for loss severity modeling found in
## Appendix A of Klugman, Panjer & Willmot (2004) and not already
## present in base R. (The generalized beta, inverse gaussian and
## log-t are not included.)
##
## In addition, the package provides "m" functions to compute
## theoretical raw moments and "lev" functions to compute limited
## moments for all the above probability laws, plus the following
## already in R: exponential, gamma, lognormal and Weibull.

## TRANSFORMED BETA FAMILY

## Tranformed beta distribution
k <- 1:10
shape1 <- 3
shape2 <- 4
shape3 <- 5
rate <- 10

showgraphs <- function(fun, par)
{
    df   <- match.fun(paste("d", fun, sep = ""))
    pf   <- match.fun(paste("p", fun, sep = ""))
    rf   <- match.fun(paste("r", fun, sep = ""))
    mf   <- match.fun(paste("m", fun, sep = ""))
    levf <- match.fun(paste("lev", fun, sep = ""))

    formals(df)[names(par)]   <- par
    formals(pf)[names(par)]   <- par
    formals(rf)[names(par)]   <- par
    formals(mf)[names(par)]   <- par
    formals(levf)[names(par)] <- par

    x <- rf(1000)
    limit <- seq(0, max(x), length = 10)
    k <- 1:10

    op <- par(mfrow = c(2, 2))

    hist(x, prob = TRUE, main = "Density")
    curve(df(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(ecdf(x), pch = "", main = "Distribution function", lwd = 2)
    curve(pf(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(k, emm(x, k), type = "l", lwd = 2, main = "Raw moments")
    lines(k, mf(k), col = "blue", lwd = 2, lty = 2)
    plot(limit, elev(x)(limit), type = "l", lwd = 2,
         main = "Limited expected value")
    lines(limit, levf(limit), col = "blue", lwd = 2, lty = 2)

    par(op)
}

showgraphs("trbeta", list(shape1 = 3, shape2 = 4, shape3 = 5, scale = 10))





par(op)
