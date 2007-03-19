### ===== actuar: an R package for Actuarial Science =====
###
### Demo of the loss distributions facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)
if(dev.cur() <= 1) get(getOption("device"))()

op <- par(ask = interactive() &&
          (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))

### A utility function to create graphs for probability laws
showgraphs <- function(fun, par, mfrow = c(2, 2))
{
    dist <- switch(fun,
                   trbeta = "TRANSFORMED BETA DISTRIBUTION",
                   genpareto = "GENERALIZED PARETO DISTRIBUTION",
                   burr = "BURR DISTRIBUTION",
                   invburr = "INVERSE BURR DISTRIBUTION",
                   pareto = "PARETO DISTRIBUTION",
                   invpareto = "INVERSE PARETO DISTRIBUTION",
                   llogis = "LOGLOGISTIC DISTRIBUTION",
                   paralogis = "PARALOGISTIC DISTRIBUTION",
                   invparalogis = "INVERSE PARALOGISTIC DISTRIBUTION",
                   trgamma = "TRANSFORMED GAMMA DISTRIBUTION",
                   invtrgamma = "INVERSE TRANSFORMED GAMMA DISTRIBUTION",
                   gamma = "GAMMA DISTRIBUTION",
                   invgamma = "INVERSE GAMMA DISTRIBUTION",
                   weibull = "WEIBULL DISTRIBUTION",
                   invweibull = "INVERSE WEIBULL DISTRIBUTION",
                   exp = "EXPONENTIAL DISTRIBUTION",
                   invexp = "INVERSE EXPONENTIAL DISTRIBUTION",
                   lnorm= "LOGNORMAL DISTRIBUTION",
                   pareto1 = "SINGLE PARAMETER PARETO DISTRIBUTION",
                   lgamma = "LOGGAMMA DISTRIBUTION")

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

    x <- rf(5000)
    limit <- seq(floor(min(x)), max(x), length = 10)
    k <- seq(1, 3, length = 6)

    op <- par(mfrow = mfrow, oma = c(0, 0, 2, 0))

    hist(x, prob = TRUE, xlim = c(0, max(x)),
         main = "Density")
    curve(df(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(ecdf(x), xlim = c(0, max(x)), pch = "",
         main = "Distribution function", lwd = 2)
    curve(pf(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
    plot(k, emm(x, k), type = "l", lwd = 2,
         main = "Raw moments")
    lines(k, mf(k), col = "blue", lwd = 2, lty = 2)
    plot(limit, elev(x)(limit), type = "l", lwd = 2,
         main = "Limited expected value")
    lines(limit, levf(limit), col = "blue", lwd = 2, lty = 2)
    title(main = dist, outer = TRUE)

    par(op)
}

###
### DATA SETS
###

## The package includes the individual dental claims and grouped
## dental claims data sets often referred to in Klugman, Panjer &
## Willmot (1998, 2004)
data(dental); dental
data(gdental); gdental


###
### PROBABILITY LAWS
###

## The package provides "d", "p", "q" and "r" functions for all the
## probability laws useful for loss severity modeling found in
## Appendix A of Klugman, Panjer & Willmot (2004) and not already
## present in base R, plus the loggamma distribution. (The generalized
## beta, inverse gaussian and log-t are not included.)
##
## In addition, the package provides "m" functions to compute
## theoretical raw moments and "lev" functions to compute limited
## moments for all the above probability laws, plus the following
## already in R: exponential, gamma, lognormal and Weibull.
##
## We illustrate the various distributions by plotting for each four
## graphs combining empirical and theoretical quantities: the pdf, the
## cdf, the first few raw moments, the limited expected value at a
## few limits.

## TRANSFORMED BETA FAMILY

## Transformed beta distribution
showgraphs("trbeta", list(shape1 = 3, shape2 = 4, shape3 = 5, scale = 10))

## Generalized Pareto distribution
showgraphs("genpareto", list(shape1 = 10, shape2 = 4, scale = 10))

## Burr distribution
showgraphs("burr", list(shape1 = 3, shape2 = 4, scale = 10))

## Inverse Burr distribution
showgraphs("invburr", list(shape1 = 3, shape2 = 6, scale = 10))

## Pareto distribution
showgraphs("pareto", list(shape = 6, scale = 10))

## Inverse Pareto distribution
showgraphs("invpareto", list(shape = 1, scale = 1))

## Loglogistic distribution
showgraphs("llogis", list(shape = 6, scale = 10))

## Paralogistic distribution
showgraphs("paralogis", list(shape = 3, scale = 10))

## Inverse paralogistic distribution
showgraphs("invparalogis", list(shape = 6, scale = 10))

## TRANSFORMED GAMMA FAMILY

## Transformed gamma distribution
showgraphs("trgamma", list(shape1 = 3, shape2 = 1, scale = 10))

## Inverse transformed gamma distribution
showgraphs("invtrgamma", list(shape1 = 3, shape2 = 2, scale = 10))

## Gamma distribution ('mgamma' and 'levgamma')
showgraphs("gamma", list(shape = 3, scale = 10))

## Inverse gamma distribution
showgraphs("invgamma", list(shape = 6, scale = 10))

## Weibull distribution ('mweibull' and 'levweibull')
showgraphs("weibull", list(shape = 1.5, scale = 10))

## Inverse Weibull distribution
showgraphs("invweibull", list(shape = 6, scale = 10))

## Exponential distribution ('mexp' and 'levexp')
showgraphs("exp", list(rate = 0.1))

## Inverse exponential distribution
showgraphs("invexp", list(rate = 1))

## OTHER DISTRIBUTIONS

## Lognormal distribution ('mlnorm' and 'levlnorm')
showgraphs("lnorm", list(meanlog = 1, sdlog = 1))

## Single parameter Pareto distribution
showgraphs("pareto1", list(shape = 5, min = 10))

## Loggamma distribution
showgraphs("lgamma", list(shapelog = 2, ratelog = 5))


###
### GROUPED DATA MANIPULATION
###

## Creation of grouped data objects
x <- grouped.data(groups = c(0, 25, 50, 100, 150, 250, 500),
                  line1 = c(30, 31, 57, 42, 65, 84),
                  line2 = c(26, 33, 31, 19, 16, 11))

## Extraction and replacement: only "[" and "[<-" are officially
## supported.
x[, 1]                                  # group boundaries
x[1]                                    # notice the difference
x[, -1]                                 # group frequencies
x[1:3,]                                 # first 3 groups
x[1, 2] <- 22; x                        # frequency replacement
x[1, 1] <- c(0, 20); x                  # boundary replacement

## A method for 'mean' exists for grouped data objects.
mean(x)

## Function 'hist' handles individual data only. We provide a method
## for grouped data.  Only the first frequencies column is considered.
hist(x[, -3])

## Function 'ogive' returns a function to compute the ogive of grouped
## data in any point, much like 'ecdf' does for individual
## data. Methods also exist to extract the group boundaries ('knots')
## and plot the ogive.
Fnt <- ogive(x)
summary(Fnt)
knots(Fnt)                              # group boundaries
Fnt(knots(Fnt))                         # ogive at group boundaries
plot(Fnt)                               # plot of the ogive


###
### EMPIRICAL MOMENTS CALCULATION
###

## Function 'emm' compute the k-th empirical moment of a sample,
## whether it is individual or grouped data.
emm(dental)                             # == mean(dental)
emm(gdental)                            # == mean(gdental)
emm(dental, order = 1:3)                # first three moments
emm(gdental, order = 1:3)               # idem

## Function 'elev' is similar to 'ecdf' and 'ogive' in that it returns
## a function to compute the empirical limited expected value (first
## limited moment) for any limit. There are methods for individual and
## grouped data.
lev <- elev(dental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function


###
### MINIMUM DISTANCE ESTIMATION
###

## Maximum likelihood estimation (for individual data) is well covered
## by 'fitdistr' in package MASS. We provide function 'mde' to fit
## models using three distance minimization techniques: Cramer-von
## Mises (for individual and grouped data), chi-square and layer
## average severity (both grouped data only). Usage (and inner
## working) is very similar to 'fitdistr'.
mde(dental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200), measure = "LAS")


###
### COVERAGE MODIFICATIONS
###

## Function 'coverage' is useful to obtain the probability density
## function (pdf) or cumulative distribution function (cdf) of a loss
## random variable under any combination of the following coverage
## modifications: ordinary or franchise deductible, policy limit,
## inflation, coinsurance. The function returned can then be used like
## any other pdf or cdf in modeling.
f <- coverage("gamma", deductible = 1, limit = 7)
curve(dgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 0.3))    # original
curve(f(x, 3, 1), xlim = c(0.01, 5.99), lty = 2, add = TRUE) # modified


par(op)
