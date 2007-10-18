### ===== actuar: an R package for Actuarial Science =====
###
### Quantiles for objects of class 'aggregateDist'
###
### AUTHORS:  Louis-Philippe Pouliot,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

quantile.aggregateDist <-
    function(x, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
             smooth = FALSE, names = TRUE, ...)
{
    label <- comment(x)

    ## The Normal and Normal Power approximations are the only
    ## continuous distributions of class 'aggregateDist'. They are
    ## therefore treated differently, using the 'base' quantile
    ## function qnorm().
    if (label == "Normal approximation")
        res <- qnorm(probs, get("mean", environment(x)),
                     sqrt(get("variance", environment(x))))
    else if (label == "Normal Power approximation")
    {
        mean <- get("mean", environment(x))
        variance <- get("variance", environment(x))
        skewness <- get("skewness", environment(x))
        ## Calling qnorm() and inverting the Normal Power 'standardization'
        res <- ifelse(probs <= 0.5, NA,
                      ((qnorm(probs) + 3/skewness)^2 - 9/(skewness^2) - 1) *
                      sqrt(variance) * skewness/6 + mean)
    }
    else
    {
        ## An empirical and discrete approach is used for
        ## 'aggregateDist' objects obtained from methods other than
        ## Normal and Normal Power.
        Fs <- get("y", environment(x))
        x <- get("x", environment(x))
        ind <- sapply(probs, function(q) match(TRUE, Fs >= q))

        res <-
            if (smooth)
            {
                h <- (Fs[ind] - probs) / (Fs[ind] - Fs[ind - 1])
                (1 - h) * x[ind - 1] + h * x[ind]
            }
            else
                x[ind]
    }

    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(res) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", wid = 1, digits = dig)
    }
    res
}
