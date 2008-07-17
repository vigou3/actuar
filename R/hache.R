### ===== actuar: an R package for Actuarial Science =====
###
### Credibility in the regression case using the Hachemeister (1975)
### model with possibly an adjustment to put the intercept at the
### barycenter of time (see Buhlmann & Gisler, 2005).
###
### AUTHORS: Xavier Milhaud, Tommy Ouellet, Vincent Goulet
### <vincent.goulet@act.ulaval.ca>

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
    ncontracts <- nrow(ratios)   	# number of contracts
    eff.ncontracts <- length(has.data)  # effective number of contracts
    p <- ncol(xreg)               	# rank (>= 2) of design matrix
    n <- NROW(xreg)                     # number of observations

    ## To put the intercept at the barycenter of time, transform the
    ## design matrix into a "weighted orthogonal" matrix.
    x <-
        if (adj.intercept)
        {
            w <- colSums(weights, na.rm = TRUE)/sum(weights.s)
            ## QR decomposition of the design matrix : X = QR
            Xwqr <- qr(xreg * sqrt(w))          # object qr
            Rw <- qr.R(Xwqr)                    # transition matrix R
            qr.Q(qr(sqrt(w) * xreg)) / sqrt(w)  ## orthogonal matrix Q
        }
        else
            xreg

    ## Fit linear model to each contract and make summary of each
    ## model for later extraction of key quantities.
    f <- function(i)
    {
    	z <-
            if (i %in% has.data)            # contract with data
                lm.wfit(x, ratios[i, ], weights[i, ])
            else                            # contract with no data
                lm.fit(x, rep.int(0, n))
        z[c("coefficients", "residuals", "weights", "rank", "qr")]
    }
    fits <- lapply(seq_len(ncontracts), f)

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

    ## "Weight" and credibility matrices are stored in p x p x
    ## ncontracts arrays. We must keep them equal to 0 (matrix)
    ## for contracts with no data. We also initialize the between
    ## variance-covariance matrix for later use.
    R <- function(z)                # from stats:::summary.lm
        chol2inv(z$qr$qr[1:p, 1:p, drop = FALSE])
    cred <- W <- array(0, c(p, p, ncontracts))
    W[, , has.data] <- sapply(fits[has.data], R)
    cred[, , has.data] <- A <- diag(p)

    ## Starting collective regression coefficients. Arithmetic
    ## average is coherent with credibility matrices as
    ## initialized above.
    coll <- rowSums(ind) / eff.ncontracts

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
    ## If the intercept is positioned at the barycenter of time,
    ## variance components are estimated just like in the
    ## Buhlmann-Straub model. This means that unbiased and iterative
    ## estimators are available.
    if (!adj.intercept)                 # intercept at time origin
    {
        ## Only iterative estimation is supported.
        if (method == "unbiased")
        {
            warning("using iterative estimators with regression at time origin")
            method <- "iterative"
        }

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
    }
    else                                # intercept at barycenter
    {
        ## We just need to compute the values on the diagonal of the
        ## variance-covariance matrix A using the Buhlmann-Straub
        ## estimators (see bstraub.R for details).
        for (i in seq_len(p))
            A[i, i] <- bvar.unbiased(ind[i, has.data], 1 / W[i, i, has.data],
                                     s2, eff.ncontracts)

        if (method == "iterative" &&
            diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
        {
            for (i in seq_len(p))
            {
                A[i, i] <-
                    if (A[i, i] > 0)
                        bvar.iterative(ind[i, has.data], 1 / W[i, i, has.data],
                                       s2, eff.ncontracts, start = A[i, i],
                                       tol = tol, maxit = maxit, echo = echo)
                    else
                        0
            }
        }

        ## Final credibility factors and estimator of the collective
        ## regression coefficients.
        for (i in seq_len(p))
        {
            if (A[i, i] > 0)
            {
                z <- cred[i, i, has.data] <- 1/(1 + (s2 * W[i, i, has.data]) / A[i, i])
                coll[i] <- drop(crossprod(z, ind[i, has.data])) / sum(z)
            }
            else
            {
                cred[i, i, ] <- numeric(ncontracts)
                w <- 1 / W[i, i, has.data]
                coll[i] <- drop(crossprod(w, ind[i, ])) / sum(w)
            }
        }
    }

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Credibility adjusted coefficients. The coefficients of the
    ## models are replaced with these values. That way, prediction
    ## will be trivial using predict.lm().
    for (i in seq_len(ncontracts))
    {
        fits[[i]]$coefficients <- coll + drop(cred[, , i] %*% (ind[, i] - coll))
        if (adj.intercept)   # change basis to come back at the basis of the beginning
            fits[[i]]$coefficients <- solve(Rw) %*% fits[[i]]$coefficients
    }
    if (adj.intercept)
    {
        ind <- solve(Rw) %*% ind    # change basis
        coll <- solve(Rw) %*% coll  # change basis
    }

    ## Results
    structure(list(means = list(coll, ind),
                   weights = list(1 / W),
                   unbiased = if (method == "unbiased") list(A, s2),
                   iterative = if (method == "iterative") list(A, s2),
                   cred = cred,
                   adj.models = fits,
                   nodes = list(ncontracts)),
              class = "hache",
              model = "regression")
}

plot.hache <- function(x, contractNo)
{
    ## Plot an object of class "hache", and particularly 3 regression lines
    ## *the collective regression line (blue) : evolution of the premium of the portfolio,
    ## *the individual one (red) : the same for the given contract,
    ## *the credibility regression line (green) : useful to check it lies between the two first
    ## I / O : object and no of the contract / the plot associated
    mu_ind <- mu_coll <- mu_cred <- numeric(ncol(weights))
    for (i in 1:ncol(weights))
    {
        mu_ind[i] <- x$means[[2]][1, contractNo] + i * x$means[[2]][2, contractNo]
        mu_coll[i] <- x$means[[1]][1] + i * x$means[[1]][2]
        mu_cred[i] <- x$adj.models[[contractNo]]$coefficients[1] +
            i * x$adj.models[[contractNo]]$coefficients[2]
    }
    ylim = range(c(mu_ind, mu_coll))
    plot(1:ncol(weights), ratios[contractNo, ], type = "p",
         ylim = extendrange(r = ylim, f = 0.05),
         xlab = "contract", ylab = "premiums")
    lines(1:ncol(weights), mu_ind, col = "red")
    lines(1:ncol(weights), mu_coll, col = "blue")
    lines(1:ncol(weights), mu_cred, col = "green")
}

predict.hache <- function(object, levels = NULL, newdata, ...)
{
    ## Prediction (credibility premiums) using predict.lm() on each of
    ## the adjusted individual models. This first requires to add a
    ## 'terms' component to each adjusted model.
    f <- function(z, ...)
    {
        z$terms <- object$terms
        predict.lm(z, ...)
    }

    sapply(object$adj.models, f, newdata = newdata)
}

print.hache <- function(x, ...)
    print.default(x)
