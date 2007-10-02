### ===== actuar: an R package for Actuarial Science =====
###
### Bühlmann-Straub credibility model calculations
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Sébastien Auclair, and Louis-Philippe Pouliot

bstraub <- function(ratios, weights, method = c("iterative", "unbiased"),
                    maxit = 100, tol = sqrt(.Machine$double.eps),
                    echo = FALSE, old.format = TRUE)
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
        stop("missing values are not in the same positions in weights and in ratios")
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
    method <- match.arg(method)

    if (method == "iterative")
    {
        if (ac > 0)
        {
            if (diff(range(weights, na.rm = TRUE)) > .Machine$double.eps^0.5)
            {
                if (echo)
                {
                    cat("Iteration\tBetween variance estimator\n")
                    exp <- expression(cat(" ", count, "\t\t ", at1 <- at,
                                          fill = TRUE))
                }
                else
                    exp <- expression(at1 <-  at)

                at <- ac
                count <- 0
                repeat
                {
                    eval(exp)

                    if (maxit < (count <- count + 1))
                    {
                        warning("maximum number of iterations reached before obtaining convergence")
                        break
                    }

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
        cred <- numeric(length(weights.s))
        ratios.zw <- ratios.ww
    }

    if (old.format)
    {
        warning("this output format is deprecated")
        structure(list(individual = ratios.w,
                       collective = ratios.zw,
                       weights = weights.s,
                       s2 = s2,
                       unbiased = ac,
                       iterative = at,
                       cred = cred),
                  class = "bstraub.old")
    }
    else
        structure(list(means = list(ratios.zw, ratios.w),
                       weights = list(if (a > 0) sum(cred) else weights.ss, weights.s),
                       unbiased = c(ac, s2),
                       iterative = if (!is.null(at)) c(at, s2),
                       cred = cred,
                       nodes = list(nrow(weights))),
                  class = "bstraub")
}

predict.bstraub.old <- function(object, levels = NULL, ...)
    object$collective + object$cred * (object$individual - object$collective)

predict.bstraub <- function(object, levels = NULL, ...)
    object$means[[1]] + object$cred * (object$means[[2]] - object$means[[1]])
