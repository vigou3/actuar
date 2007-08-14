### ===== actuar: an R package for Actuarial Science =====
###
### Credibility in the Regression Case
###
### The Hachemeister Regression Model (1975).
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

hache <- function(X, Y, weights, TOL = 1E-6, echo = FALSE)
{
    Call <- match.call()

    ## X = data matrix (I x T)
    ## Y = design matrix (T x N)
    I <- nrow(X)     # number of contracts
    T <- ncol(X)     # number of years
    N <- ncol(Y)     # dimension of design matrix
    
    ## Coefficients for each contract
    beta <- matrix(0, nrow = I, ncol = N)
    for (k in 1:I) beta[k, ] = coef(lm(X[k, ] ~ Y[, -1], weights = weights[k, ]))
    
    ## Estimation of s^2.
    s2 <- 0
    for (k in 1:I)
    {
        diff <- X[k, ] - Y %*% beta[k, ]
        s2 <- s2 + t(diff) %*% diag(weights[k, ]) %*% diff
    }
    s2 <- as.vector( s2 / (I * (T - N)) )

    ## Diagonal matrices filled with diagonal elements one are used
    ## as starting values for the Zi matrices. The matrix of each
    ## contract being of length N x N, the resulting matrices are
    ## stocked in an array of length N x N x I where I is the
    ## number of contracts.
    Z <- array(as.vector(diag(N)), dim = c(N, N, I))
    
    ## First estimation of betaTotal.
    t1 = t2 = 0
    for (k in 1:I)
    {
        t1 = t1 + Z[, , k]
        t2 = t2 + Z[, , k] %*% beta[k, ]
    }
    betaTotal <- solve(t1) %*% t2
    
    ## Iterative estimation of the parameters.
    repeat
    {
        ## Transferring the value of betaTotal in another variable.
        ## This value will be used as the stop-criterion.
        if (echo) print(as.vector(betaTotal))
        betaTotal2 <- betaTotal 

        ## Estimation of A.
        A <- 0
        for (k in 1:I)
        {
            diff <- beta[k, ] - betaTotal
            A <- A + Z[, , k] %*% diff %*% t(diff)
        }
        A <- A / (I - 1)

        ## Matrix A being symmetrical, A is replaced by ( A + t(A) ) / 2.
        A <-  ( A + t(A) ) / 2

        ## Estimation of the Zi matrices.
        for (k in 1:I) Z[, , k] <- A %*% solve( ( A + s2 * solve( t(Y) %*% diag(weights[k, ]) %*% Y ) ) )

        ## New estimation of betaTotal.
        t1 = t2 = 0
        for (k in 1:I)
        {
            t1 = t1 + Z[, , k]
            t2 = t2 + Z[, , k] %*% beta[k, ]
        }
        betaTotal <- solve(t1) %*% t2
        
        ## Stop-criterion.
        if (max(abs( (betaTotal - betaTotal2) / betaTotal2 ) ) < TOL) break
    }

    ## Final estimation of A and Z with the final betaTotal.
    A <- 0
    for (k in 1:I)
    {
        diff <- beta[k, ] - betaTotal
        A <- A + Z[, , k] %*% diff %*% t(diff)
    }
    A <- A / (I - 1)
    A <-  ( A + t(A) ) / 2

    for (k in 1:I) Z[, , k] <- A %*% solve( ( A + s2 * solve( t(Y) %*% diag(weights[k, ]) %*% Y ) ) )

    ## Credibility-adjusted estimator for beta.
    betaAdj <- matrix(0, nrow = I, ncol = N)
    for (k in 1:I) betaAdj[k, ] = Z[, , k] %*% beta[k, ] + (diag(N) - Z[, , k]) %*% betaTotal

    ## Credibility premiums
    p <- matrix(0, nrow = I, ncol = T)
    for (k in 1:I) p[k, ] = t( Y %*% betaAdj[k, ] )

    ## Results
    structure(list(beta = beta,
                   betaAdj = betaAdj,
                   betaTotal = betaTotal,
                   cred = Z,
                   s2 = s2,
                   a = A,
                   p = p,
                   call = Call),
              class = "hache")
}

print.hache <- function(x, ...)
{
    I <- nrow(x$beta)
    res <- matrix(NA, I + 1, 4)
    res[, 1:2] <- rbind(x$beta, t(x$betaTotal))
    res[, 3:4] <- rbind(x$betaAdj, NA)
    colnames(res) <- c("Intercept", "Slope", "Adj. intercept", "Adj. slope")
    rownames(res) <- c(1:I, "Total")

    cat("\nCall: ", deparse(x$call), "\n\n")

    print(res)

    cat("\nWithin contract variance: ", x$s2, "\n\n")
    cat("Between contract variance:\n")
    print(x$a)
}
    
