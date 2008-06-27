### ===== actuar: an R package for Actuarial Science =====
###
### Credibility in the regression case using the Hachemeister (1975)
### model with possibly an adjustment to put the intercept at the
### barycenter of time (see Buhlmann & Gisler, 2005).
###
### AUTHORS: Xavier Milhaud, Vincent Goulet <vincent.goulet@act.ulaval.ca>

hache <- function(ratios, weights, xreg, adj.intercept = FALSE,
                  method = c("unbiased", "iterative"),
                  tol = sqrt(.Machine$double.eps),
                  maxit = 100, echo = FALSE)
{
    Call <- match.call()
    method <- match.arg(method)

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

    ## Frequently used values. Note that 'xreg' is guaranteed to have
    ## at least 2 columns when called from cm.
    weights.s <- rowSums(weights, na.rm = TRUE) # contract total weights
    has.data <- which(weights.s > 0)	# contracts with data
    all.na <- which(weights.s == 0)     # contracts without data
    ncontracts <- nrow(ratios)   	# number of contracts
    eff.ncontracts <- sum(weights.s > 0) # effective number of contracts
    p <- ncol(xreg)               	# rank (>= 2) of design matrix
    n <- NROW(xreg)                     # number of observations
    s <- seq_len(ncontracts)	        # sequence 1, ..., ncontracts
    #k <- seq_len(ncol(weights))	        # sequence 1, ..., n

    ## If there should be an adjustment to take the barycenter of time
    ## as intercept in the regression, it is computed as the *across
    ## contract* weighted average of the regressor. Only one regressor
    ## (i.e. simple regression) is supported in this case.
    if (adj.intercept)
    {
        if (p > 2 && adj.intercept)
            stop("only one regressor is supported in regression at the barycenter of time")
        j0 <- sum(weights %*% xreg[, 2], na.rm = TRUE)/sum(weights.s)
        xreg[, 2] <- xreg[, 2] - j0
    }

    ## Fit linear model to each contract and make summary of each
    ## model for later extraction of key quantities.
    f <- function(i)
    {
    	if (i %in% all.na)              # contract with no data
            lm.fit(xreg, rep.int(0, n))
        else                            # contract with data
            lm.wfit(xreg, ratios[i, ], weights[i, ])
    }
    fits <- lapply(s, f)

    ## Individual regression coefficients
    ind <- sapply(fits, coef)

    ## Individual variance estimators
    S <- function(z)                    # from stats:::summary.lm
    {
        r <- z$residuals
        w <- z$weights
        if (is.null(w))
            rss <- sum(r^2)
        else
            rss <- sum(w * r^2)
        rss/(n - p)
    }
    sigma2 <- sapply(fits, S)

    ## === ESTIMATION OF WITHIN VARIANCE ===
    s2 <- sum(sigma2) / eff.ncontracts

    ## === ESTIMATION OF THE BETWEEN VARIANCE-COVARIANCE MATRIX ===
    ##
    ## The procedure is different depending if there is an adjustment
    ## of the intercept or not.
    ##
    ## If there is no adjustment, we use an iterative procedure
    ## similar to the Bischel-Straub estimator. Following Goovaerts &
    ## Hoogstad (1987), the stopping criterion is based on the
    ## collective regression coefficients estimates.
    ##
    ## If the intercept is positioned at the barycenter of time, then
    ## intercept and slope variance components are estimated just like
    ## in the Buhlmann-Straub model. This means that unbiased and
    ## iterative estimators are available.
    if (!adj.intercept)                 # intercept at time origin
    {
        ## Only iterative estimation is supported.
        if (method == "unbiased")
            warning("using iterative estimators with regression at time origin")

        ## "Weight" and credibility matrices are stored in p x p x
        ## ncontracts arrays. We must keep them equal to 0 (matrix)
        ## for contracts with no data. We also initialize the between
        ## variance-covariance matrix A for later use.
        R <- function(z)                # from stats:::summary.lm
            chol2inv(z$qr$qr[1:p, 1:p, drop = FALSE])
        cred <- W <- array(0, c(p, p, ncontracts))
        W[, , has.data] <- sapply(fits[has.data], R)
        cred[, , has.data] <- A <- diag(p)

        ## Starting collective regression coefficients. Arithmetic
        ## average is coherent with above credibility matrices.
        coll <- rowSums(ind) / eff.ncontracts

        ## If printing of iterations was asked for, start by printing a
        ## header and the starting values.
        if (echo)
        {
            cat("Iteration\tCollective regression coefficients\n")
            exp <- expression(cat(" ", count, "\t\t ", coll1 <- coll,
                fill = TRUE))
        }
        else
            exp <- expression(coll1 <-  coll)

        ## Iterative procedure
        count <- 0
        repeat
        {
            eval(exp)

            ## Stop after 'maxit' iterations
            if (maxit < (count <- count + 1))
            {
                warning("maximum number of iterations reached before obtaining convergence")
                break
            }

            ## Calculation of the between variance-covariance matrix.
            A[] <- rowSums(sapply(has.data, function(i)
                                cred[, , i] %*% tcrossprod(ind[, i] - coll))) /
                                    (eff.ncontracts - 1)

            ## Symmetrize A
            A <- (A + t(A))/2

            ## New credibility matrices
            cred[, , has.data] <- sapply(has.data, function(i)
                                         A %*% solve(A + s2 * W[, , i]))

            ## New collective regression coefficients
            cred.s <- apply(cred[, , has.data], c(1, 2), sum)
            coll <- solve(cred.s,
                          rowSums(sapply(has.data, function(i)
                                         cred[, , i] %*% ind[, i])))

            ## Test for convergence
            if (max(abs((coll - coll1)/coll1)) < tol)
                break
        }

        ## Final calculation of the between variance-covariance matrix and
        ## credibility matrices.
        A[] <- rowSums(sapply(has.data, function(i)
                              cred[, , i] %*% tcrossprod(ind[, i] - coll))) /
                                  (eff.ncontracts - 1)
        A <- (A + t(A))/2
        cred[, , has.data] <- sapply(has.data, function(i)
                                     A %*% solve(A + s2 * W[, , i]))

        ## Credibility adjusted coefficients. The coefficients of the
        ## models are replaced with these values. That way, prediction
        ## will be trivial using predict.lm().
        for (i in s)
            fits[[i]]$coefficients <-
                coll + drop(cred[, , i] %*% (ind[, i] - coll))
    }
    else                                # intercept at barycenter
    {
        g <- function(i)
        {
            if (!all(is.na(weights[i,])))
                (1 / weights.s[i]) * sum(k * weights[i,], na.rm = TRUE)
            else
                0
        }
        ## compute vector of individual centres of gravity
        barycentres <- sapply(s, g)
        ## compute all d_i
        d <- numeric(nrow(ratios))
        for (i in s)
        {
            if (!all(is.na(weights[i,])))
            {
                for (j in k)
                    d[i] <- d[i] + (weights[i,j] / weights.s[i]) * (j-barycentres[i])^2
            }
        }
        ## compute adjusted weights w*
        weightsadj <- d*weights.s

        tau0c <- bvar.unbiased(ind[1,], weights.s, s2, ncontracts)
        tau0c <- max(tau0c, 0)
        tau1c <- bvar.unbiased(ind[2,], weightsadj, s2, ncontracts)
        tau1c <- max(tau1c, 0)

        if (method == "iterative")
        {
            tau0t <-
                if (tau0c > 0)
                {
                    if (diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
                        bvar.iterative(ind[1,], weights.s, s2, ncontracts, start = tau0c,
                                       tol = tol, maxit = maxit, echo = echo)
                    else
                        tau0c
                }
                else
                    0
            tau0 <- tau0t

            tau1t <-
                if (tau1c > 0)
                {
                    if (diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
                        bvar.iterative(ind[2,], weightsadj, s2, ncontracts, start = tau1c,
                                       tol = tol, maxit = maxit, echo = echo)
                    else
                        tau1c
                }
                else
                    0
            tau1 <- tau1t
        }
        else
        {
            tau0 <- tau0c
            tau1 <- tau1c
            tau0t <- NULL
            tau1t <- NULL
        }

        ## final calculation of param a11, a22, beta0 and beta1
        a11 <- weights.s / (weights.s + s2/tau0)
        if (tau0 == 0)
            beta0coll <- sum(weights.s * ind[1,]) / sum(weights.s)
        else
            beta0coll <- sum(a11 * ind[1,]) / sum(a11)

        a22 <- weights.s / (weights.s + s2/tau1)
        if (tau1 == 0)
            beta1coll <- sum(weightsadj * ind[2,]) / sum(weightsadj)
        else
            beta1coll <- sum(a22 * ind[2,]) / sum(a22)

        ## to have the same format of results
        coll <- c(beta0coll, beta1coll)
        ## fill in the credibility matrices with a11 and a22 coef.
        for (i in s)
            cred[, , i] <- diag(c(a11[i], a22[i]), p, p)
        ## Credibility estimators of beta_0 and beta_1 : adjusted coefficients
        for (i in s)
        {
            fits[[i]]$coefficients[1] <- coll[1] + a11[i] * (ind[1,i] - coll[1])
            fits[[i]]$coefficients[2] <- coll[2] + a22[i] * (ind[2,i] - coll[2])
        }
    }

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Results
    structure(list(means = list(coll, ind),
                   weights = if (!adjust) list(cred.s, lapply(covUnscaled[has.data], solve)) else NULL,
                   unbiased = if (!adjust) NULL else { if (method == "unbiased") list(tau0, tau1) else NULL },
                   iterative = if (!adjust) list(A, s2) else { if (method == "iterative") list(tau0, tau1) },
                   cred = cred,
                   adj.models = fits,
                   nodes = list(ncontracts)),
              class = "hache",
              model = "regression")
}

predict.hache <- function(object, levels = NULL, newdata, ...)
{
    ## Predictors can be given as a simple vector for one dimensional
    ## models. For use in predict.lm(), these must be converted into a
    ## data frame.
    if (is.null(dim(newdata)))
        newdata <- data.frame(xreg = newdata)

    ## Prediction (credibility premiums) using predict.lm() on each of
    ## the adjusted individual models.
    sapply(object$adj.models, predict, newdata = newdata)
}

print.hache <- function(x, ...)
    print.default(x)
