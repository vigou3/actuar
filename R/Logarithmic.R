### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,r}logarithmic functions to compute
### characteristics of the logarithmic distribution. The version used
### in these functions has probability mass function
###
###   Pr[X = x] = -p^x/(x log(1 - p)), x = 1, 2, ...
###
### This is the standard parametrization in the literature; see for
### example https://en.wikipedia.org/wiki/Logarithmic_distribution.
###
### NOTE: Klugman, Panjer & Willmot (Loss Models) introduce the
### logarithmic distribution as a limiting case of the zero truncated
### negative binomial. In this setting, parameter 'p' above would be
### the probability of *failure* (a.k.a q) of the zero truncted
### negative binomial.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dlogarithmic <- function(x, prob, log = FALSE)
    .External("actuar_do_dpq", "dlogarithmic", x, prob, log)

plogarithmic <- function(q, prob, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "plogarithmic", q, prob, lower.tail, log.p)

## qlogarithmic <- function(p, prob, lower.tail = TRUE, log.p = FALSE)
##     .External("actuar_do_dpq", "qlogarithmic", p, prob, lower.tail, log.p)

rlogarithmic <- function(n, prob)
    .External("actuar_do_random", "rlogarithmic", n, prob)

rlogarithmicLS <- function(n, prob)
{
    p <- 1 - prob
    t <- -p/log(prob)
    x <- 1
    u <- runif(1)
    while (u > t)
    {
        u <- u - t
        x <- x + 1
        t <- t * p * (x - 1)/x
    }
    x
}

rlogarithmicLK <- function(n, prob)
{
    r <- log(prob)
    x <- 1
    v <- runif(1)

    if (v >= 1 - prob) return(x)

    u <- runif(1)
    q <- 1 - exp(r * u)
    if (v <= q^2) return(round(1 + log(v)/log(q)))
    if (v <= q)   return(1)
    2
}
