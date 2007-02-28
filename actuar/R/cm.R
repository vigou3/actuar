### ===== actuar: an R package for Actuarial Science =====
###
### Credibility Models
###
### Fit a credibility model in the formulation of variance components
### as described in Dannenburg, Kaas and Goovaerts (1996). Models
### supported are part of a generalized hierarchical credibility
### theory as introduced in Dannenburg (1995).
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Louis-Philippe Pouliot

cm <- function(formula, data, years, weights, subset, TOL = 1E-6, echo = FALSE)
{
    cl <- match.call()

    ## Check if 'formula' expresses hierarchical interactions.
    if (any(duplicated(attr(terms(formula), "order"))))
        stop("unsupported interactions in 'formula'")

    ## A subset of 'data' is created, and it is then divided into
    ## three data frames: one for the portfolio structure, one for the ratios
    ## and one for the weights.

    levs <- rev(rownames(attr(terms(formula), "factors")))

    if (missing(subset))
        r <- TRUE
    else
    {
        e <- substitute(subset)
        r <- eval(e, data, parent.frame())
        if (!is.logical(r))
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
        #if (length(unique(data[[struct]])) < length(r[r]))    ####### probleme... a voir
         #   stop("Hierarchical conflict in 'formula'/'subset'")
    }

    if (!inherits(data, "data.frame"))
        stop("'data' must be a 'data frame'")

    if (missing(years))
    {
        if (!missing(weights))
            stop("ratios have to be specified if weigths are")
        years <- suppressWarnings(names(data)[!(names(data) %in% levs)])
    }
    else
    {
        nl <- as.list(1:ncol(data))
        names(nl) <- names(data)
        years <- eval(substitute(years), nl, parent.frame())
    }


    ## Check if interactions are consistent with the data

    nstruct <- c(length(years), sapply(levs, function(x) length(unique(data[[x]]))), pf = 1)
    if (!all(sort(nstruct, decreasing = TRUE) == nstruct))
        stop("hierarchical interactions are inconsistent with the data")

    ratios <- data[r, years]

    ## If weights are not specified, use equal weights.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing ratios not allowed when the matrix of weights is not specified")
        wt <- matrix(1, nrow(data[r, ]), length(years))
    }
    else
    {
        if (ncol(ratios) < 2)
            stop("there must be at least one contract with at least two years of experience")
        if (nrow(data) < 2)
            stop("there must be more than one contract")

        weights <- eval(substitute(weights), nl, parent.frame())
        if (length(weights) != length(years))
            stop("'years' and 'weights' must have the same length")
        wt <- data[r, weights]
        if (!identical(which(is.na(ratios)), which(is.na(wt))))
            stop("missing values are inconsistent in 'ratios'/'weights' data")
    }

    ## Coerce the structural part of the data frame to class 'factor', with
    ## levels considering the order of the first occurence of a value in the frame.

    ## Bind a column of '1's representing the affiliation to the one global portfolio.
    ## Used to symmetrize further calculations.

    data <- cbind(pf = 1, data[r, ])
    levs <- c(levs, "pf")
    nLevels <- length(levs)
    fstruct <- data.frame(sapply(data[levs], function(x) factor(x, levels = unique(x))))


    ## A list expressing the affiliation structure
    aff <- vector("list", nLevels - 1)
    for (i in 1:(nLevels - 1))
        aff[[i]] <- tapply(fstruct[[levs[i + 1]]], fstruct[[levs[i]]], function(x) unique(x))


    ind.weight <- rowSums(wt, na.rm = TRUE)
    ind.means <- rowSums(ratios * wt, na.rm = TRUE) / ind.weight

    ## s^2
    s2 <- sum(wt * (ratios - ind.means) ^ 2, na.rm = TRUE)

    denoms <- numeric(nLevels)
    for (i in 2:nLevels) denoms[i] <- nstruct[i] - nstruct[i + 1]

    s2 <- s2 /(denoms[1] <- sum(!is.na(ratios)) - nstruct[[2]])

    ## Create vectors for values to be outputted.

    param <- rep(s2, nLevels)        # The structure parameters
    cred <- vector("list", nLevels - 1)      # The credibility factors
    w. <- vector("list", nLevels)    # The credibility weights
    M <- vector("list", nLevels)     # The individual and collective estimators

    ## Avoid evaluating argument 'echo' at every iteration below
    if (echo)
        exp <- expression(print(paramt <- param))
    else
        exp <- expression({paramt <- param})

    ## Iterative estimation of the structure parameters
    repeat
    {
        eval(exp)

        ## Individual estimators are initialized at every iteration.
        weight <- ind.weight
        means <- ind.means

        for (i in 1:(nLevels - 1))
        {
            cred[[i]] <- 1/(1 + param[i]/(param[i + 1] * weight))
            weight. <- tapply(cred[[i]], aff[[i]], sum)
            means. <- tapply(cred[[i]] * means, aff[[i]], sum) / weight.
            param[i + 1] <- sum(cred[[i]] *
                            (means - rep(means., table(aff[[i]]))[order(aff[[i]])])^2) / denoms[i + 1]

            w.[[i]] <- weight
            weight <- weight.
            M[[i]] <- means
            means <- means.
        }
        p <- ifelse(any(param <= TOL), which(param > TOL), TRUE)
        if (max(abs((param[p] - paramt[p])/paramt[p])) < TOL)
                break
    }
    w.[[nLevels]] <- weight. ## is it necessary?
    M[[nLevels]] <- means.
    res <- list(param = param,
                weights = w., ## is it?
                means = M,
                cred = cred,
                call = cl,
                data = data,
                levs = levs,
                aff = aff)
    class(res) <- "cm"
    res
}

print.cm <- function(x, ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium: ", x$means[[length(x$means)]], "\n")
    cat("  Expected variance: ", x$param[1],"\n")
    for (i in (l <- 2:length(x$param)))
        cat(" ", letters[i-1], " :", x$param[i], "\n")
    for (i in (l-1))
    {
        a <- sapply(x$aff[[i]], function(z) as.character(unique(x$data[x$levs][[i + 1]]))[z])
        cat("\nLevel: ", x$levs[i], "\n")
        m <- data.frame(x$mean[[i]], x$cred[[i]], a)
        dimnames(m) <- list(unique(unlist(x$data[[x$levs[i]]])),
                            c("ind. estimator", "cred. factor", paste(x$levs[i + 1], "affiliation")))
        print(m)
    }
}
