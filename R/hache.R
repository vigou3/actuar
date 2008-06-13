### ===== Actuar: an R package for Actuarial Science =====
###
### Credibility in the Regression Case
###
### The Hachemeister Regression Model (1975).
###
### AUTHORS: Xavier Milhaud, Vincent Goulet <vincent.goulet@act.ulaval.ca>

hache <- function(ratios, weights, xreg, method = c("unbiased", "iterative"),
                       adjust = FALSE, tol = sqrt(.Machine$double.eps),
                       maxit = 100, echo = FALSE)
{
    Call <- match.call()
    method <- match.arg(method)

    ## If weights are not specified, use equal weights as in
    ## Buhlmann's model.
    if (missing(weights))
    {
      	stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
    }

    ## Check other bad arguments.
    if (ncol(ratios) < 2)
        stop("there must be at least one node with more than one period of experience")
    if (nrow(ratios) < 2)
        stop("there must be more than one node")
    if (!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in weights and in ratios")
    if (all(!weights, na.rm = TRUE))
        stop("no available data to fit model")
    if (method == "unbiased" && !adjust)
        warning("estimators are calculated iteratively, argument method should be iterative")

    ## Frequently used values
    weights.s <- rowSums(weights, na.rm = TRUE)
    has.data <- which(weights.s > 0)	# contracts with data
    all.na <- which(weights.s == 0)     # contracts without data
    ncontracts <- sum(weights.s > 0)   	# effective number of contracts
    p <- NCOL(xreg) + 1               	# dimension of design matrix
    s <- seq_len(nrow(ratios))	        # vector 1:5 with Hachemeister data, number of contracts
    k <- seq_len(ncol(weights))	        # vector 1:12 with Hachemeister data, periods

    if (p > 2 && adjust)
        stop("if adjust = TRUE, regressor must be vector and not a multidimensionnal object")

    xreg <- cbind(xreg)                 # force dims and colnames
    ## Compute centre of gravity of time j0 in data
    if (adjust)
    {
        j0 <- drop(crossprod(colSums(weights, na.rm = TRUE), k))/sum(weights.s)
        xreg <- xreg - j0
    }

    ## Fit linear model to each contract and make summary of each
    ## model for later extraction of key quantities.
    fo <- as.formula(paste("y ~ ", paste(colnames(xreg), collapse = "+")))
    f <- function(i)
    {
    	if (all(is.na(ratios[i,])))
    	{
            ratios[i,] <- 0
            DF <- data.frame(y = ratios[i, ], xreg, w = weights[i, ])
            lm(fo, data = DF, weights = NULL)
    	}
        else
        {
            DF <- data.frame(y = ratios[i, ], xreg, w = weights[i, ])
            lm(fo, data = DF, weights = w)
        }
    }
    fits <- lapply(s, f)
    sfits <- lapply(fits, summary)

    ## Regression coefficients, residuals and the analogue of the inverse
    ## of the total contract weights (to be used to compute the
    ## credibility matrices). for each contract
    ind <- sapply(fits, coef)
    sigma2 <- sapply(sfits, "[[", "sigma")^2

    ## === ESTIMATION OF WITHIN VARIANCE ===
    s2 <- sum(sigma2) / ncontracts

    ## Starting credibility matrices, stored in an array
    ## of dimension p x p x ncontracts.
    cred <- array(diag(p), dim = c(p, p, nrow(ratios))) # identity matrices
    cred[, , all.na] <- matrix(0, p, p)  # coherent with above cred. matrices

    if (!adjust)
    {
        ## covUnscaled equivalent to alpha_j factor in credibility factor ("Surveys of Actuarial studies",p.55)
        covUnscaled <- vector("list", nrow(ratios))
        covUnscaled[has.data] <- lapply(sfits[has.data],"[[", "cov.unscaled")
        covUnscaled[all.na] <- list(matrix(0, nrow = p, ncol = p))

        ## === ESTIMATION OF THE BETWEEN VARIANCE-COVARIANCE MATRIX ===
        ##
        ## This is an iterative procedure similar to the Bischel-Straub
        ## estimator. Following Goovaerts & Hoogstad, stopping criterion
        ## is based in the collective regression coefficients estimates.
        ##
        ## Starting collective regression coefficients.
        coll <- rowSums(ind) / ncontracts

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
            A <- rowSums(sapply(has.data, function(i) cred[, , i] %*% tcrossprod(ind[, i] - coll))) / (ncontracts - 1)
            dim(A) <- c(p, p)

            ## Symmetrize A
            A <- (A + t(A))/2

            ## New credibility matrices
            cred[, , has.data] <- sapply(covUnscaled[has.data], function(w) A %*% solve(A + s2 * w))
            dim(cred) <- c(p, p, nrow(ratios))

            ## New collective regression coefficients
            cred.s <- apply(cred[, , has.data], c(1, 2), sum)
            coll <- solve(cred.s,
                          rowSums(sapply(has.data, function(i) cred[, , i] %*% ind[, i])))

            ## Test for convergence
            if (max(abs((coll - coll1)/coll1)) < tol)
                break
        }

        ## Final calculation of the between variance-covariance matrix and
        ## credibility matrices.
        A <- rowSums(sapply(has.data, function(i) cred[, , i] %*% tcrossprod(ind[, i] - coll))) / (ncontracts - 1)
        dim(A) <- c(p, p)
        A <- (A + t(A))/2
        cred[, , has.data] <- sapply(covUnscaled[has.data], function(w) A %*% solve(A + s2 * w))
        dim(cred) <- c(p, p, nrow(ratios))

        ## Credibility adjusted coefficients. The coefficients of the
        ## models are replaced with these values. That way, prediction
        ## will be trivial using predict.lm().
        for (i in s)
            fits[[i]]$coefficients <- coll + drop(cred[, , i] %*% (ind[, i] - coll))
    }
    else
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
