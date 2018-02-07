### ===== actuar: An R Package for Actuarial Science =====
###
### Pure bayesian credibility calculations.
###
### AUTHORS: Alexandre Parent <alexandre.parent.12@ulaval.ca>,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

bayes <- function(x, likelihood =
                         c("poisson", "bernoulli", "geometric",
                           "exponential", "normal", "binomial",
                           "negative binomial", "gamma",
                           "pareto"),
                  shape, rate = 1, scale = 1/rate,
                  shape1, shape2, mean = 0, sd = 1,
                  size, shape.lik, sd.lik, min)
{
    likelihood <- match.arg(likelihood)

    ## We need to treat separately the (Single Parameter, or
    ## Translated) Pareto/Gamma case given the different form of the
    ## individual mean and the "credibility factor" (which isn't one,
    ## really).
    if (likelihood == "pareto")
    {
        if (missing(min))
            stop("lower bound of the likelihood missing")
        if (missing(shape) || (missing(rate) && missing(scale)))
            stop("one of the Gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
        coll <- shape * scale
        vars <- c(NA, NA)                # not pertinent here

        ## Computation of individual means and credibility factors
        ## differs depending on the type of data provided in argument.
        if (is.null(x))                 # no data
            cred <- ind.means <- n <- 0
        else if (is.vector(x, mode = "numeric")) # atomic vector
        {
            n <- length(x)
            sumlog <- sum(log(x)) - n * log(min)
            ind.means <- n/sumlog
            cred <- 1/(1 + 1/(scale * sumlog))
        }
        else                            # matrix or data frame
        {
            n <- ncol(x)
            sumlog <- rowSums(log(x)) - n * log(min)
            ind.means <- n/sumlog
            cred <- 1/(1 + 1/(scale * sumlog))
        }
    }
    ## Now the usual linear Bayes cases.
    else
    {
        if (likelihood == "poisson")
        {
            if (missing(shape) || (missing(rate) && missing(scale)))
                stop("one of the Gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
            coll <- shape * scale
            vars <- c(coll * scale, coll)
            K <- 1/scale
        }
        else if (likelihood == "bernoulli")
        {
            if (missing(shape1) || missing(shape2))
                stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
            K <- shape1 + shape2
            coll <- shape1/K
            vars <- (shape1 * shape2) * c(1, K)/(K^2 * (K + 1))
        }
        else if (likelihood == "binomial")
        {
            if (missing(shape1) || missing(shape2))
                stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
            if (missing(size))
                stop("parameter \"size\" of the likelihood missing")
            K <- (shape1 + shape2)/size
            coll <- shape1/K
            vars <- (shape1 * shape2) * c(1, K)/(K^2 * (shape1 + shape2 + 1))
        }
        else if (likelihood == "geometric")
        {
            if (missing(shape1) || missing(shape2))
                stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
            K <- shape1 - 1
            coll <- shape2/K
            vars <- shape2 * (shape1 + shape2 - 1)/(K * (K - 1))
            vars <- c(vars/K, vars)
        }
        else if (likelihood == "negative binomial")
        {
            if (missing(shape1) || missing(shape2))
                stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
            if (missing(size))
                stop("parameter \"size\" of the likelihood missing")
            K <- (shape1 - 1)/size
            coll <- shape2/K
            vars <- shape2 * (shape1 + shape2 - 1)/(K * (shape1 - 2))
            vars <- c(vars/K, vars)
        }
        else if (likelihood == "exponential")
        {
            if (missing(shape) || (missing(rate) && missing(scale)))
                stop("one of the Gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
            K <- shape - 1
            coll <- 1/(K * scale)
            vars <- c(coll^2, coll/scale)/(shape - 2)
        }
        else if (likelihood == "gamma")
        {
            if (missing(shape) || (missing(rate) && missing(scale)))
                stop("one of the gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
            if (missing(size))
                stop("parameter \"shape.lik\" of the likelihood missing")
            K <- (shape - 1)/shape.lik
            coll <- 1/(K * scale)
            vars <- c(coll^2, coll/scale)/(shape - 2)
        }
        else if (likelihood == "normal")
        {
            if (missing(sd.lik))
                stop("standard deviation of the likelihood missing")
            coll <- mean
            vars <- c(sd, sd.lik)^2
            K <- vars[2L]/vars[1L]
        }
        else
            stop("unsupported likelihood")

        ## Computation of individual means and credibility factors
        ## differs depending on the type of data provided in argument.
        if (is.null(x))                 # no data
            cred <- ind.means <- n <- 0
        else if (is.vector(x, mode = "numeric")) # atomic vector
        {
            n <- length(x)
            ind.means <- mean(x)
            cred <- n/(n + K)
        }
        else                            # matrix or data frame
        {
            n <- ncol(x)
            ind.means <- rowMeans(x)
            cred <- n/(n + K)
        }
    }

    structure(list(means = list(coll, ind.means),
                   weights = rep_len(1, n),
                   unbiased = vars,
                   iterative = NULL,
                   cred = cred,
                   nodes = 1L),
              class = "bayes",
              model = "Linear Bayes")

}

## Premium calculation is identical to the Buhlmann-Straub case; no
## need for another method. See bstraub.R for the definition.
# predict.bayes <- predict.bstraub
