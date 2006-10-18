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

### A utily function to create graphs for probability laws
showgraphs <- function(fun, par, mfrow = c(2, 2))
{
    dist <- switch(fun,
                   trbeta = "TRANSFORMED BETA DISTRIBUTION",
                   burr = "BURR DISTRIBUTION",
                   llogis = "LOGLOGISTIC DISTRIBUTION",
                   paralogis = "PARALOGISTIC DISTRIBUTION",
                   genpareto = "GENERALIZED PARETO DISTRIBUTION",
                   pareto = "PARETO DISTRIBUTION",
                   pareto1 = "SINGLE PARAMETER PARETO DISTRIBUTION",
                   invburr = "INVERSE BURR DISTRIBUTION",
                   invpareto = "INVERSE PARETO DISTRIBUTION",
                   invparalogis = "INVERSE PARALOGISTIC DISTRIBUTION",
                   trgamma = "TRANSFORMED GAMMA DISTRIBUTION",
                   invtrgamma = "INVERSE TRANSFORMED GAMMA DISTRIBUTION",
                   invgamma = "INVERSE GAMMA DISTRIBUTION",
                   invweibull = "INVERSE WEIBULL DISTRIBUTION",
                   invexp = "INVERSE EXPONENTIAL DISTRIBUTION",
                   lgamma = "LOGGAMMA DISTRIBUTION",
                   lnorm= "LOGNORMAL DISTRIBUTION",
                   gamma = "GAMMA DISTRIBUTION",
                   exp = "EXPONENTIAL DISTRIBUTION",
                   weibull = "WEIBULL DISTRIBUTION")

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
    k <- seq(1, 5, length = 10)

    op <- par(mfrow = mfrow, oma = c(0, 0, 2, 0))

    hist(x, prob = TRUE, main = "Density")
    curve(df(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(ecdf(x), pch = "", main = "Distribution function", lwd = 2)
    curve(pf(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(k, emm(x, k), type = "l", lwd = 2, main = "Raw moments")
    lines(k, mf(k), col = "blue", lwd = 2, lty = 2)
    plot(limit, elev(x)(limit), type = "l", lwd = 2,
         main = "Limited expected value")
    lines(limit, levf(limit), col = "blue", lwd = 2, lty = 2)
    title(main = dist, outer = TRUE)

    par(op)
}


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

## Transformed beta distribution
showgraphs("trbeta", list(shape1 = 3, shape2 = 4, shape3 = 5, scale = 10))

## Generalized Pareto distribution
showgraphs("genpareto", list(shape1 = 3, shape2 = 4, scale = 10))

## Burr distribution
showgraphs("burr", list(shape1 = 3, shape2 = 4, scale = 10))

## Inverse Burr distribution
showgraphs("invburr", list(shape1 = 3, shape2 = 4, scale = 10))

## Pareto distribution
showgraphs("pareto", list(shape = 3, scale = 10))

## Inverse Pareto distribution
showgraphs("invpareto", list(shape = 0.01, scale = 10))

## Paralogistic distribution
showgraphs("paralogis", list(shape = 3, scale = 10))

## Inverse paralogistic distribution
showgraphs("invparalogis", list(shape = 3, scale = 10))

## Loglogistic distribution
showgraphs("llogis", list(shape = 3, scale = 10))


## TRANSFORMED GAMMA FAMILY

## Transformed gamma distribution
showgraphs("trgamma", list(shape1 = 3, shape2 = 4, scale = 10))

## Gamma distribution ('mgamma' and 'levgamma')
showgraphs("gamma", list(shape = 3, scale = 10))

## Weibull distribution ('mweibull' and 'levweibull')
showgraphs("weibull", list(shape = 3, scale = 10))

## Exponential distribution ('mexp' and 'levexp')
showgraphs("exp", list(rate = 0.1))


## INVERSE TRANSFORMED GAMMA FAMILY

## Inverse transformed gamma distribution
showgraphs("invtrgamma", list(shape1 = 3, shape2 = 4, scale = 10))

## Inverse gamma distribution
showgraphs("invgamma", list(shape = 3, scale = 10))

## Inverse Weibull distribution
showgraphs("invweibull", list(shape = 3, scale = 10))

## Inverse exponential distribution
showgraphs("invexp", list(rate = 0.1))


## OTHER DISTRIBUTIONS

## Lognormal distribution ('mlnorm' and 'levlnorm')
showgraphs("lnorm", list(meanlog = 5, sdlog = 1))

## Single parameter Pareto distribution
showgraphs("pareto1", list(shape = 3, min = 10))


###
### GROUPED DATA MANIPULATION
###

## Creation of grouped data objects
x <- grouped.data(group = c(0, 25, 50, 100, 150, 250, 500),
                  line1 = c(30, 31, 57, 42, 65, 84),
                  line2 = c(26, 33, 31, 19, 16, 11))

## Extraction and replacement: only "[" and "[<-" are officially
## supported
x[, 1]                                  # group boundaries
x[1]                                    # notice the difference
x[, -1]                                 # group frequencies
x[1:3,]                                 # first 3 groups
x[1, 2] <- 22; x                        # frequency replacement
x[1, 1] <- c(0, 20); x                  # boundary replacement


par(op)
