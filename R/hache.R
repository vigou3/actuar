### ===== actuar: an R package for Actuarial Science =====
###
### Credibility in the regression case using the Hachemeister (1975)
### model with possibly an adjustment to put the intercept at the
### barycenter of time (see Buhlmann & Gisler, 2005).
###
### AUTHORS: Xavier Milhaud, Tommy Ouellet, Vincent Goulet
### <vincent.goulet@act.ulaval.ca>

hache <- function(ratios, weights, formula, data, adj.intercept = FALSE,
                  method = c("unbiased", "iterative"),
                  tol = sqrt(.Machine$double.eps),
                  maxit = 100, echo = FALSE)
{
    Call <- match.call()

    ## If weights are not specified, use equal weights as in
    ## Buhlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
    }

    ## Check other bad arguments.
    if (NCOL(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")
    if (NROW(ratios) < 2)
        stop("there must be more than one node")
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in 'weights' and in 'ratios'")
    if (all(!weights, na.rm = TRUE))
        stop("no available data to fit model")

    ## Build the design matrix
    mf <- model.frame(formula, data, drop.unused.levels = TRUE)
    mt <- attr(mf, "terms")
    xreg <- model.matrix(mt, mf)

    ## Do computations in auxiliary functions.
    res <-
        if (adj.intercept)
            hache.barycenter(ratios, weights, xreg,
                             method = match.arg(method),
                             tol = tol, maxit = maxit, echo = echo)
        else
            hache.origin(ratios, weights, xreg,
                         tol = tol, maxit = maxit, echo = echo)

    ## Add the terms object to the result for use in predict.hache()
    ## [and thus predict.lm()].
    res$terms <- mt

    ## Results
    attr(res, "class") <- "hache"
    attr(res, "model") <- "regression"
    res
}

predict.hache <- function(object, levels = NULL, newdata, ...)
{
    ## If model was fitted at the barycenter of time (there is a
    ## transition matrix in the object), then also convert the
    ## regression coefficients in the base of the (original) design
    ## matrix.
    if (!is.null(R <- object$transition))
    {
        for (i in seq_along(object$adj.models))
        {
            b <- coefficients(object$adj.models[[i]])
            object$adj.models[[i]]$coefficients <- solve(R, b)
        }
    }

    ## Prediction (credibility premiums) using predict.lm() on each of
    ## the adjusted individual models. This first requires to add a
    ## 'terms' component to each adjusted model.
    f <- function(z, ...)
    {
        z$terms <- object$terms
        unname(predict.lm(z, ...))
    }

    sapply(object$adj.models, f, newdata = newdata)
}

print.hache <- function(x, ...)
    print.default(x)

## This function does not plot the prevision because it uses the
## curve() function, and we'd just want 1 prediction at t = 13
## with Hachemeister data, which is impossible to obtain
## because of arguments 'from', 'to' and 'n' and the fact that
## curve() uses a fonction of x and not a vector as first parameter.
plot.hache <- function(x, contractNo, from = NULL, to = NULL, n = 101,
                       type = c("predictions", "heterogeneity"), add = FALSE,
                       main = NULL)
{
    ## check the class of the object
    if (!inherits(x, "hache"))
                stop("use only with \"hache\" objects")

    type <- match.arg(type)
    if (type == "predictions")
    {
        ## Plot object of class "hache", more particularly 3 lines:
        ## * the collective regression line (blue)
        ## * the individual one (red)
        ## * the credibility regression line (green)

        ## If model was fitted at the barycenter of time (there is a
        ## transition matrix in the object), then also convert the
        ## regression coefficients in the base of the (original) design
        ## matrix.
        if (!is.null(R <- x$transition))
        {
            x$means[[1]] <- solve(R, x$means[[1]])
            x$means[[2]] <- solve(R, x$means[[2]])
            for (i in seq_along(x$adj.models))
            {
                b <- coefficients(x$adj.models[[i]])
                x$adj.models[[i]]$coefficients <- solve(R, b)
            }
        }

        mu_ind <- mu_coll <- mu_cred <- numeric(n)
        ## degree of the regression (ex : quadratic => degree = 2)
        degree <- length(x$means[[1]]) - 1
        abscisses <- seq(from = from, to = to, length.out = n)
        comp.premium.individual <- function(abscisses)
        {
            ## compute values of the individual regression line
            for (i in 1:n)
            {
                mu_ind[i] <- x$means[[2]][1, contractNo]
                for (j in 1:degree)
                    mu_ind[i] <- mu_ind[i] + abscisses[i]^j * x$means[[2]][j+1, contractNo]
            }
            mu_ind
        }
        comp.premium.collective <- function(abscisses)
        {
            for (i in 1:n)
            {
                mu_coll[i] <- x$means[[1]][1]
                for (j in 1:degree)
                    mu_coll[i] <- mu_coll[i] + abscisses[i]^j * x$means[[1]][j+1]
            }
            mu_coll
        }
        comp.premium.credibility <- function(abscisses)
        {
            for (i in 1:n)
            {
                mu_cred[i] <- x$adj.models[[contractNo]]$coefficients[1]
                for (j in 1:degree)
                    mu_cred[i] <- mu_cred[i] + abscisses[i]^j * x$adj.models[[contractNo]]$coefficients[j+1]
            }
            mu_cred
        }
        mu_ind <- comp.premium.individual(abscisses)
        mu_coll <- comp.premium.collective(abscisses)
        ylim = range(c(mu_ind, mu_coll))
        ## curve() requires a function of 'abscisses' in first param
        curve(comp.premium.individual, from = from, to = to, n = n, col = "red",
              add = add, xlab = "time", ylab = "premiums", main = main, ylim = ylim)
        curve(comp.premium.collective, from = from, to = to, n = n, col = "blue",
              add = TRUE, xlab = "time", ylab = "premiums", main = main, ylim = ylim)
        curve(comp.premium.credibility, from = from, to = to, n = n, col = "green",
              add = TRUE, xlab = "time", ylab = "premiums", main = main, ylim = ylim)
        legend("topright", legend = c("collective premium", "individual premium", "credibility premium"),
               text.width = strwidth("credibility premium"), lty = c(1,1,1), col = c("blue","red","green"),
               xjust = 1, yjust = 1, cex = 0.6, title = "Regression lines")
    }
    else
    {
        ## construct the structure of the summary, generated in our fashion.
        ## -> 'stats' argument : changes from classical boxplot().
        ## Meaning of each column :
        ##  * 1st and 2nd element : ind. mean - 0.5*sqrt(within variance of the contract)
        ##    (2 equal elements so as to hide the whiskers)
        ##  * 3rd element : ind. mean of the contract
        ##  * 4th and 5th element : ind. mean + 0.5*sqrt(within variance of the contract)
        statsIntercept <- matrix(0, nrow = 5, ncol = nrow(x$ratios))
        statsSlope <- matrix(0, nrow = 5, ncol = nrow(x$ratios))
        for (j in 1:nrow(x$ratios))
        {
            statsIntercept[1, j] <- statsIntercept[2, j] <- x$means[[2]][1,j] - 0.5*sqrt(x$within[j])
            statsSlope[1, j] <- statsSlope[2, j] <- x$means[[2]][2,j] - 0.5*sqrt(x$within[j])
            statsIntercept[3, j] <- x$means[[2]][1,j]
            statsSlope[3, j] <- x$means[[2]][2,j]
            statsIntercept[4, j] <- statsIntercept[5, j] <- x$means[[2]][1,j] + 0.5*sqrt(x$within[j])
            statsSlope[4, j] <- statsSlope[5, j] <- x$means[[2]][2,j] + 0.5*sqrt(x$within[j])
        }
        ## we have to specify different attributes to bxp()
        ## 'n' attribute = number of observations
        n <- numeric(nrow(x$ratios))
        for (i in 1:nrow(x$ratios))
            n[i] <- length(x$ratios[i, ])
        ## 'confidence interval' attribute
        conf <- matrix(NA, nrow = 2, ncol = nrow(x$ratios))
        out <- numeric()          ## 'out' attribute
        group <- numeric()        ## 'group' attribute
        ## 'names' attribute
        names <- vector("character", nrow(x$ratios))
        for (i in 1:nrow(x$ratios))
            names[i] <- i

        bx.pInter <- list(stats = statsIntercept, n = n, conf = conf, out = out,
                     group = group, names = names)
        bx.pSlope <- list(stats = statsSlope, n = n, conf = conf, out = out,
                     group = group, names = names)
        ## graphs : first the graph of homogeneity in intercepts
        par(mfrow = c(1,2))
        bxp(bx.pInter, show.names = TRUE, medcol = "red", medlwd = 7,
            boxwex = 0.01, staplewex = 30, ylab = "Within variance",
            main = "Homogeneity of the portfolio \n (intercept)", xlab = "contract No")
        arrows(x0 = 1.5, y0 = min(statsIntercept[3, ]), x1 = 1.5, y1 = max(statsIntercept[3, ]),
               length = 0.1, angle = 15, code = 3, lty = 7)
        text(1.5, 0.5*(min(statsIntercept[3, ])+max(statsIntercept[3, ])), "between variance",
             pos = 4, cex = 0.7, srt = 90)
        legend("topright", legend = "individual premium",
               text.width = strwidth("individual premium"), lty = -1, pch = 19, col = "red",
               xjust = 1, yjust = 1, cex = 0.7)
        ## then the graph of homogeneity in slopes
        bxp(bx.pSlope, show.names = TRUE, medcol = "red", medlwd = 7,
            boxwex = 0.01, staplewex = 30, ylab = "Within variance",
            main = "Homogeneity of the portfolio \n (slope)", xlab = "contract No")
        arrows(x0 = 1.5, y0 = min(statsSlope[3, ]), x1 = 1.5, y1 = max(statsSlope[3, ]),
               length = 0.1, angle = 15, code = 3, lty = 7)
        text(1.5, 0.5*(min(statsSlope[3, ])+max(statsSlope[3, ])), "between variance",
             pos = 4, cex = 0.7, srt = 90)
        legend("topright", legend = "individual premium",
               text.width = strwidth("individual premium"), lty = -1, pch = 19, col = "red",
               xjust = 1, yjust = 1, cex = 0.7)
    }
}
