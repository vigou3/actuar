### ===== actuar: An R Package for Actuarial Science =====
###
### Pure bayesian credibility calculations.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

bayesian <- function(x, likelihood =
                            c("poisson", "bernoulli", "geometric",
                              "exponential", "normal",
                              "binomial", "negative binomial", "gamma"),
                     shape, rate = 1, scale = 1/rate,
                     shape1, shape2,
                     sd.lik = 1, mean = 0, sd = 1)
{
    likelihood <- match.arg(likelihood)

    if (!is.vector(x, "numeric"))
        stop ("data must be a numeric vector for Bayesian models")

    if (likelihood == "bernoulli")
    {
        if (missing(shape1) || missing(shape2))
            stop("one of the Beta prior parameter \"shape1\" or \"scale2\" missing")
        coll = shape1/(shape1 + shape2)
        K = shape1 + shape2
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
        vars = c(sd, sd.lik)
        K = sd.lik/sd
    }
    else
        stop("unsupported likelihood")

    n <- length(x)
    structure(list(means = c(coll, if (length(x) > 0) mean(x) else 0),
                   weights = rep_len(1, n),
                   unbiased = vars,
                   iterative = NULL,
                   cred = n/(n + K),
                   nodes = 1L),
              class = "bayesian",
              model = "Pure Bayesian")

}

## Premium calculation is identical to the Buhlmann-Straub case; no
## need for another method.
predict.bayesian <- predict.bstraub
