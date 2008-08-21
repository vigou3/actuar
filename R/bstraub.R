### ===== actuar: an R package for Actuarial Science =====
###
### Bühlmann-Straub credibility model calculations.
###
### Computation of the between variance estimators has been moved to
### external functions bvar.unbiased() and bvar.iterative() to share
### with hache().
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair, Louis-Philippe Pouliot

bstraub <- function(ratios, weights, method = c("unbiased", "iterative"),
                    tol = sqrt(.Machine$double.eps), maxit = 100,
                    echo = FALSE, old.format = FALSE)
{
    ## If weights are not specified, use equal weights as in
    ## Bühlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
    }

    ## Check other bad arguments.
    if (ncol(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")
    if (nrow(ratios) < 2)
        stop("there must be more than one node")
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in 'weights' and in 'ratios'")
    if (all(!weights, na.rm = TRUE))
        stop("no available data to fit model")

    ## Individual weighted averages. It could happen that a contract
    ## has no observations, for example when applying the model on
    ## claim amounts. In such a situation, we will put the total
    ## weight of the contract and the weighted average both equal to
    ## zero. That way, the premium will be equal to the credibility
    ## weighted average, as it should, but the contract will have no
    ## contribution in the calculations.
    weights.s <- rowSums(weights, na.rm = TRUE)
    ratios.w <- ifelse(weights.s > 0, rowSums(weights * ratios, na.rm = TRUE) / weights.s, 0)

    ## Size of the portfolio.
    nobs <- ncol(ratios)
    ncontracts <- sum(weights.s > 0)
    ntotal <- sum(!is.na(weights))

    ## Collective weighted average.
    weights.ss <- sum(weights.s)

    ## Estimators of individual within variances
    sigma2 <- rowSums(weights * (ratios - ratios.w)^2) / (nobs-1)

    ## Estimation of s^2 : mean of within variances.
    s2 <-  sum(weights * (ratios - ratios.w)^2, na.rm = TRUE) / (ntotal - ncontracts)

    ## First estimation of a. Always compute the unbiased estimator.
    a <- bvar.unbiased(ratios.w, weights.s, s2, ncontracts)

    ## Iterative estimation of a. Compute only if
    ## 1. asked to in argument;
    ## 2. weights are not all equal (Bühlmann model).
    ## 3. the unbiased estimator is > 0;
    method <- match.arg(method)

    if (method == "iterative" &&
        diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
    {
        a <-
            if (a > 0)
                bvar.iterative(ratios.w, weights.s, s2, ncontracts, start = a,
                               tol = tol, maxit = maxit, echo = echo)
            else
                0
    }

    ## Final credibility factors and estimator of the collective mean.
    if (a > 0)
    {
        cred <- 1 / (1 + s2/(weights.s * a))
        ratios.zw <- drop(crossprod(cred, ratios.w)) / sum(cred)
    }
    else
    {
        cred <- numeric(length(weights.s))
        ratios.zw <- drop(crossprod(weights.s, ratios.w)) / sum(weights.s)
    }

    if (old.format)
    {
        warning("this output format is deprecated")
        structure(list(individual = ratios.w,
                       collective = ratios.zw,
                       weights = weights.s,
                       s2 = s2,
                       unbiased = if (method == "unbiased") a,
                       iterative = if (method == "iterative") a,
                       cred = cred),
                  class = "bstraub.old",
                  model = "Buhlmann-Straub")
    }
    else
        structure(list(ratios = ratios,
                       means = list(ratios.zw, ratios.w),
                       weights = list(if (a > 0) sum(cred) else weights.ss, weights.s),
                       unbiased = if (method == "unbiased") c(a, s2),
                       iterative = if (method == "iterative") c(a, s2),
                       within = sigma2,
                       cred = cred,
                       nodes = list(nrow(weights))),
                  class = "bstraub",
                  model = "Buhlmann-Straub")
}

predict.bstraub.old <- function(object, ...)
    object$collective + object$cred * (object$individual - object$collective)

predict.bstraub <- function(object, levels = NULL, newdata, ...)
    object$means[[1]] + object$cred * (object$means[[2]] - object$means[[1]])

plot.bstraub <- function(x, contractNo, add = FALSE, main = NULL,
                         type = c("predictions","heterogeneity"))
{
    ## check the class of the object
    if (!inherits(x, "bstraub"))
        stop("use only with \"bstraub\" objects")

    type <- match.arg(type)
    if (type == "predictions")
    {
        ## plot an object of class 'bstraub' and more particularly :
        ##   * the collective regression line (blue) : constant,
        ##   * the individual regression line (red) : constant,
        ##   * the prediction of the credibility premium : constant.

        ## draw the data required in the model
        plot(1:(ncol(x$ratios)+1), c(x$ratios[contractNo, ],NA), type = "p",
             xlab = "time", ylab = "premiums",
             main = paste("Evolution of the premiums: contract ", contractNo,"(B-S model)"))
        ## add the collective mean
        lines(1:ncol(x$ratios), rep(x$means[[1]],ncol(x$ratios)),
              type = "l", col = "blue")
        ## add the individual mean
        lines(1:ncol(x$ratios), rep(x$means[[2]][contractNo],ncol(x$ratios)),
              type = "l", col = "red")
        ## add prediction to the graph : "star point"
        points(ncol(x$ratios)+1,
               predict(x, newdata = data.frame(ncol(x$ratios)+1))[contractNo],
               pch = 8, col = "green")
        legend("topright", legend = c("observations", "collective premium", "individual premium", "credibility premium prediction"),
               text.width = strwidth("credibility premium prediction"), lty = c(-1,1,1,-1),
               pch = c(21,-1,-1,8), col = c("black","blue","red","green"), xjust = 1, yjust = 1, cex = 0.6)
    }
    else
    {
        ## construct the structure of the summary, generated in our fashion.
        ## -> 'stats' argument : changes from classical boxplot()
        ## each column has the following meaning
        ##  * 1st and 2nd element : ind. mean - 0.5*sqrt(within variance of the contract)
        ##    (2 equal elements so as to hide the whiskers)
        ##  * 3rd element : ind. mean of the contract
        ##  * 4th and 5th element : ind. mean + 0.5*sqrt(within variance of the contract)
        stats <- matrix(0, nrow = 5, ncol = nrow(x$ratios))
        for (j in 1:nrow(x$ratios))
        {
            stats[1, j] <- stats[2, j] <- x$means[[2]][j] - 0.5*sqrt(x$within[j])
            stats[3, j] <- x$means[[2]][j]
            stats[4, j] <- stats[5, j] <- x$means[[2]][j] + 0.5*sqrt(x$within[j])
        }
        ## specify different attributes necessary to bxp()
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

        bx.p <- list(stats = stats, n = n, conf = conf, out = out,
                     group = group, names = names)
        bxp(bx.p, show.names = TRUE, medcol = "red", medlwd = 7,
            boxwex = 0.01, staplewex = 30, ylab = "Within variance",
            main = "Homogeneity of the portfolio", xlab = "contracts ")
        arrows(x0 = 1.5, y0 = min(stats[3, ]), x1 = 1.5, y1 = max(stats[3, ]),
               length = 0.1, angle = 15, code = 3, lty = 2)
        text(1.5, 0.5*(min(stats[3, ])+max(stats[3, ])), "between variance", pos = 4,
             cex = 0.7, srt = 90)
        legend("topright", legend = "individual premium",
               text.width = strwidth("individual premium"), lty = -1, pch = 19,
               col = "red", xjust = 1, yjust = 1, cex = 0.7)
    }
}

bvar.unbiased <- function(x, w, within, n)
{
    w.s <- sum(w)
    x.w <- drop(crossprod(w, x)) / w.s

    w.s * (drop(crossprod(w, (x - x.w)^2)) - (n - 1) * within) / (w.s^2 - sum(w^2))
}

bvar.iterative <- function(x, w, within, n, start,
                           tol = sqrt(.Machine$double.eps), maxit = 100,
                           echo = FALSE)
{
    if (echo)
    {
        cat("Iteration\tBetween variance estimator\n")
        exp <- expression(cat(" ", count, "\t\t ", a1 <- a, fill = TRUE))
    }
    else
        exp <- expression(a1 <-  a)

    a <- start
    count <- 0

    repeat
    {
        eval(exp)

        if (maxit < (count <- count + 1))
        {
            warning("maximum number of iterations reached before obtaining convergence")
            break
        }

        cred <- 1 / (1 + within/(w * a))
        x.z <- drop(crossprod(cred, x)) / sum(cred)
        a <- drop(crossprod(cred, (x - x.z)^2)) / (n - 1)

        if (abs((a - a1)/a1) < tol)
            break
    }
    a
}
