### ===== actuar: An R Package for Actuarial Science =====
###
### Pure bayesian credibility calculations.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

bayes <- function(x, likelihood =
                         c("poisson", "bernoulli", "geometric",
                           "exponential", "normal",
                           "binomial", "negative binomial", "gamma"),
                  shape, rate = 1, scale = 1/rate,
                  shape1, shape2,
                  sd.lik = 1, mean = 0, sd = 1)
{
    likelihood <- match.arg(likelihood)

    if (likelihood == "bernoulli")
    {
        if (missing(shape1) || missing(shape2))
            stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
        K = shape1 + shape2
        coll = shape1/K
        vars = (shape1 * shape2) * c(1, K)/(K^2 * (K + 1))
    }
    else if (likelihood == "poisson")
    {
        if (missing(shape) || (missing(rate) && missing(scale)))
            stop("one of the gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
        coll = shape * scale
        vars = c(coll * scale, coll)
        K = 1/scale
    }
    else if (likelihood == "geometric")
    {
        if (missing(shape1) || missing(shape2))
            stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
        a <- shape1
        b <- shape2
        K <- a - 1
        coll = b/K
        vars <- b * (a + b - 1)/(K * (K - 1))
        vars <- c(vars/K, vars)
    }
    else if (likelihood == "exponential")
    {
        if (missing(shape) || (missing(rate) && missing(scale)))
            stop("one of the gamma prior parameter \"shape\", \"rate\" or \"scale\" missing")
        K <- shape - 1
        coll = 1/(K * scale)
        vars = c(coll/scale, coll^2)/(shape - 2)
    }
    else if (likelihood == "normal")
    {
        if (missing(sd.lik))
            stop ("standard deviation of the likelihood missing")
        coll = mean
        vars = c(sd, sd.lik)^2
        K = var[2L]/vars[1L]
    }
    else
        stop("unsupported likelihood")

    ## In the pure Bayesian case we allow NULL data (makes sense
    ## since there is no estimation to be carried), an atomic vector
    ## (for a single contract), or a matrix or data frame (in which
    ## case we compute the Bayesian premium for each contract).
    ## Computation of the number of years 'n' and the individual means
    ## differ if data is a vector or a matrix/data frame.
    if (is.null(dim(x)))
    {
        n <- length(x)
        ind.means <- if (n > 0) mean(x) else 0
    }
    else
    {
        n <- ncol(x)
        ind.means <- rowMeans(x)
    }

    structure(list(means = list(coll, ind.means),
                   weights = rep_len(1, n),
                   unbiased = vars,
                   iterative = NULL,
                   cred = n/(n + K),
                   nodes = 1L),
              class = "bayes",
              model = "Pure Bayesian")

}

## Premium calculation is identical to the Buhlmann-Straub case; no
## need for another method. See bstraub.R for the definition.
# predict.bayes <- predict.bstraub
