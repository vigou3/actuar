### ===== Actuar: an R package for Actuarial Science =====
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

    ## Frequently used values
    weights.s <- rowSums(weights, na.rm = TRUE)
    has.data <- which(weights.s > 0)	# contracts with data
    all.na <- which(weights.s == 0)     # contracts without data
    ncontracts <- sum(weights.s > 0)   	# effective number of contracts
    p <- NCOL(xreg) + 1               	# dimension of design matrix
    s <- seq_len(nrow(ratios))	        # vector 1:5 with Hachemeister datas, real number of contracts
    k <- seq_len(ncol(weights))	        # vector 1:12 with Hachemeister datas, periods

    ## Fit linear model to each contract and make summary of each
    ## model for later extraction of key quantities.
    xreg <- cbind(xreg)                 # force dims and colnames
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
    ## Starting credibility matrices and collective regression
    ## coefficients. The credibility matrices are stored in an array
    ## of dimension p x p x ncontracts.
    cred <- array(diag(p), dim = c(p, p, nrow(ratios))) # identity matrices
    cred[, , all.na] <- matrix(0, p, p)
    coll <- rowSums(ind) / ncontracts         				# coherent with above cred. matrices

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

    ## Add names to the collective coefficients vector.
    names(coll) <- rownames(ind)

    ## Results
    structure(list(means = list(coll, ind),
                   weights = list(cred.s, lapply(covUnscaled[has.data], solve)),
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
