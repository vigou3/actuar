### ===== actuar: an R package for Actuarial Science =====
###
### Quantiles for objects of class 'aggregateDist'
###
### AUTHORS:  Louis-Philippe Pouliot, Vincent Goulet <vincent.goulet@act.ulaval.ca>

quantile.aggregateDist <-
    function(x, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
             approx.lin = FALSE, names = TRUE, ...)
{
    label <- comment(x)

    ## The Normal and Normal Power approximations are the only
    ## continuous distributions of class 'aggregateDist'. They are
    ## therefore treated differently, using the 'base' quantile
    ## function qnorm().
    if (label == "Normal approximation")
        ## A call to qnorm() with the given the moments of the distribution
        qnorm(probs, get("mean" ,environment(x)), sqrt(get("var", environment(x))))
    if (label == "Normal Power approximation")
    {
        mean <- get("mean" ,environment(x))
        var <- get("var", environment(x))
        skewness <- get("skewness", environment(x))
        ## Calling qnorm() and inverting the Normal Power 'standardization'
        ans <- ifelse(probs <= 0.5, NA,
                      ((qnorm(probs) + 3/skewness)^2
                       - 9/(skewness^2) - 1) * sqrt(var) * skewness/6 + mean)
    }
    else
    {
        ## An empirical and discrete approach is used for 'aggregateDist'
        ## objects issued of methods other than Normal or Normal Power.
        Fs <- get("y", environment(x))
        x <- get("x", environment(x))
        upper <- sapply(probs, function(q) x[min(which(Fs > q))])

        ans <- upper
        if (approx.lin)
        {
            lower <- sapply(probs, function(q) max(which(Fs < q)))
            h <- (Fs[upper] - probs) / (Fs[upper] - Fs[lower])
            ans <- (1 - h) * lower + h * upper
        }
    }
    if (names)
    {
        dig <- max(2, getOption("digits"))
        names(ans) <- formatC(paste(100 * probs, "%", sep = ""),
                              format = "fg", wid = 1, digits = dig)
    }
    ans
}
