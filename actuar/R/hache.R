### ===== actuar: an R package for Actuarial Science =====
###
### Credibility in the Regression Case
###
### The Hachemeister Regression Model (1975).
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

hache <- function(ratios, weights, xreg, tol = sqrt(.Machine$double.eps),
                  maxit = 100, echo = FALSE)
{
    Call <- match.call()

    ## Frequently used values
    ncontracts <- NROW(ratios)          # number of contracts
    nyears <- NCOL(ratios)              # number of years
    p <- NCOL(xreg) + 1                 # dimension of design matrix
    s <- seq_len(ncontracts)            # 1:ncontracts

    ## Fit linear model to each contract and make summary of each
    ## model for later extraction of key quantities.
    f <- function(i)
    {
        y <- ratios[i, ]
        w <- weights[i, ]
        lm(y ~ xreg, weights = w)
    }
    fits <- lapply(s, f)
    sfits <- lapply(fits, summary)

    ## Regression coefficients, residuals and the analog of the inverse
    ## of the total contract weights (to be used to compute the
    ## credibility matrices). for each contract
    ind <- sapply(fits, coef)
    r <- sapply(fits, residuals)
    sigma2 <- sapply(sfits, "[[", "sigma")^2
    weights.s <- lapply(sfits, "[[", "cov.unscaled")

    ## === ESTIMATION OF WITHIN VARIANCE ===
    s2 <- mean(sigma2)

    ## === ESTIMATION OF THE BETWEEN VARIANCE-COVARIANCE MATRIX ===
    ##
    ## This is an iterative procedure similar to the Bischel-Straub
    ## estimator. Following Goovaerts & Hoogstad, stopping criterion
    ## is based in the collective regression coefficients estimates.
    ##
    ## Starting credibility matrices and collective regression
    ## coefficients. The credibility matrices are stored in an array
    ## of dimension p x p x ncontracts.
    cred <- array(diag(p), dim = c(p, p, ncontracts)) # identity matrices
    coll <- rowMeans(ind)         # coherent with above cred. matrices

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

        ## As calculated here, the between variance-covariance matrix
        ## is actually a vector. It is turned into a matrix by adding
        ## a 'dim' attribute.
        A <- rowSums(sapply(s,
                            function(i) cred[, , i] %*% tcrossprod(ind[, i] - coll))) / (ncontracts - 1)
        dim(A) <- c(p, p)

        ## Symmetrize A
        A <- (A + t(A))/2

        ## New credibility matrices
        cred <- sapply(weights.s, function(w) A %*% solve(A + s2 * w))
        dim(cred) <- c(p, p, ncontracts)

        ## New collective regression coefficients
        cred.s <- apply(cred, c(1, 2), sum)
        coll <- solve(cred.s,
                      rowSums(sapply(s, function(i) cred[, , i] %*% ind[, i])))

        ## Test for convergence
        if (max(abs((coll - coll1)/coll1)) < tol)
            break
    }

    ## Final calculation of the between variance-covariance matrix and
    ## credibility matrices.
    A <- rowSums(sapply(s,
                        function(i) cred[, , i] %*% tcrossprod(ind[, i] - coll))) / (ncontracts - 1)
    dim(A) <- c(p, p)
    A <- (A + t(A))/2
    cred <- sapply(weights.s, function(w) A %*% solve(A + s2 * w))
    dim(cred) <- c(p, p, ncontracts)

    ## Credibility adjusted coefficients. The coefficients of the
    ## models are replaced with these values. That way, prediction
    ## will be trivial using predict.lm().
    for (i in s)
        fits[[i]]$coefficients <- coll + drop(cred[, , i] %*% (ind[, i] - coll))

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Results
    structure(list(means = list(coll, ind),
                   weights = list(cred.s, lapply(weights.s, solve)),
                   unbiased = NULL,
                   iterative = list(A, s2),
                   cred = cred,
                   adj.models = fits,
                   nodes = list(ncontracts)),
              class = "hache",
              model = "regression")
}

predict.hache <- function(object, levels = NULL, newdata, ...)
{
    ## Predictors can be given as a simple vector for one dimensional
    ## models. For use in predict.lm(), these must be tranformed into
    ## a data frame.
    if (is.null(dim(newdata)))
        newdata <- data.frame(xreg = newdata)

    ## Prediction (credibility premiums) using predict.lm() on each of
    ## the adjusted individual models.
    sapply(object$adj.models, predict, newdata = newdata)
}

print.hache <- function(x, ...)
    print.default(x)

## print.hache <- function(x, ...)
## {
##     ncontracts <- nrow(x$beta)
##     res <- matrix(NA, ncontracts + 1, 4)
##     res[, 1:2] <- rbind(x$beta, t(x$betaTotal))
##     res[, 3:4] <- rbind(x$betaAdj, NA)
##     colnames(res) <- c(" Intercept", " Slope", " Adj. intercept", " Adj. slope")
##     rownames(res) <- c(rownames(res, do.NULL = FALSE, prefix = "Contract.")[-(ncontracts + 1)], "Total")

##     cat("\nCall: ", deparse(x$call), "\n\n")

##     print(res)

##     cat("\nWithin contract variance: ", x$s2, "\n\n")
##     cat("Between contract variance:\n")
##     print(x$a)
## }

##     rownames(premiums) <- rownames(premiums, do.NULL = FALSE, prefix = "Contract.")
##     colnames(premiums) <- colnames(premiums, do.NULL = FALSE, prefix = "Time.")
##     print(premiums)
