### ===== actuar: an R package for Actuarial Science =====
###
### Bühlmann-Straub credibility model calculations
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair, and Louis-Philippe Pouliot

bstraub <- function(ratios, weights,
                    heterogeneity = c("iterative", "unbiased"),
                    tol = 1E-6, echo = FALSE)
{
    Call <- match.call()

    ## If weights are not specified, use equal weights as in
    ## Bühlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when weights are not supplied")
        weights <- array(1, dim(ratios))
        model <- "Buhlmann"
    }
    else
        model <- "Buhlmann-Straub"

    ## Check other bad arguments.
    if (ncol(ratios) < 2)
        stop("there must be at least one contract with more than one period of experience")
    if (nrow(ratios) < 2)
        stop("there must be more than one contract")
    if(!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in weights and in ratios")

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
    ncontracts <- sum(weights.s > 0)
    ntotal <- sum(!is.na(weights))

    ## Collective weighted average.
    weights.ss <- sum(weights.s)
    ratios.ww <- sum(weights.s * ratios.w) / weights.ss

    ## Estimation of s^2.
    s2 <-  sum(weights * (ratios - ratios.w)^2, na.rm = TRUE) / (ntotal - ncontracts)

    ## First estimation of a. Always compute the unbiased estimator.
    ac <- weights.ss * (sum(weights.s * (ratios.w - ratios.ww)^2) - (ncontracts - 1) * s2) / (weights.ss^2 - sum(weights.s^2))

    ## Iterative estimation of a. Compute only if
    ## 1. asked to in argument;
    ## 2. the unbiased estimator is > 0;
    ## 3. weights are not all equal (Bühlmann model).
    heterogeneity <- match.arg(heterogeneity)

    if (heterogeneity == "iterative")
    {
        if (ac > 0)
        {
            if (diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
            {
                if (echo)
                    exp <- expression(print(at1 <-  at))
                else
                    exp <- expression(at1 <-  at)

                at <- ac
                repeat
                {
                    eval(exp)

                    cred <- 1 / (1 + s2/(weights.s * at))
                    ratios.zw <- sum(cred * ratios.w) / sum(cred)
                    at <- sum(cred * (ratios.w - ratios.zw)^2) / (ncontracts - 1)

                    if (abs((at - at1)/at1) < tol)
                        break
                }
            }
            else
                at <- ac
        }
        else
            at <- 0
        a <- at
    }
    else
    {
        a <- ac
        at <- NULL
    }

    ## Final credibility factors and estimator of the collective mean.
    if (a > 0)
    {
        cred <- 1 / (1 + s2/(weights.s * a))
        ratios.zw <- sum(cred * ratios.w) / sum(cred)
    }
    else
    {
        cred <- 0
        ratios.zw <- ratios.ww
    }

    structure(list(model = model,
                   individual = ratios.w,
                   collective = ratios.zw,
                   weights = weights.s,
                   ncontracts = ncontracts,
                   s2 = s2,
                   cred = cred,
                   call = Call,
                   unbiased = ac,
                   iterative = at),
              class = "bstraub")
}

print.bstraub <- function(x, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:       ", x$collective, "\n")
    cat("  Within contract variance: ", x$s2, "\n")
    cat("  Portfolio heterogeneity:  ",
        if (is.null(x$iterative)) x$unbiased else x$iterative, "\n\n")
    invisible(x)
}

predict.bstraub <- function(object, ...)
    object$collective + object$cred * (object$individual - object$collective)

summary.bstraub <- function(object, ...)
    structure(object, class = c("summary.bstraub", class(object)), ...)

print.summary.bstraub <- function(x, ...)
{
    cat("\nCredibility model:", x$model, "\n\n")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:        ", x$collective, "\n")
    cat("  Within contract variance:  ", x$s2, "\n")
    cat("  Portfolio heterogeneity:   ", x$unbiased, " (unbiased)\n")
    cat("                             ", if (!is.null(x$iterative)) x$iterative else NA, " (iterative)\n")
    cat("  Credibility constant:      ",
        x$s2 / if (is.null(x$iterative)) x$unbiased else x$iterative,
        "\n\n")
    cat("Detailed premiums\n\n")
    cred <- cbind(1:x$ncontracts, x$individual, x$weights,
                  x$cred, predict(x))
    colnames(cred) <- c(" Contract", "Ind. premium", "Weight",
                        "Cred. factor", "Cred. premium")
    rownames(cred) <- rep("", x$ncontracts)
    print(cred, ...)
    invisible(x)
}
