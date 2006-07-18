cm <- function(formula, data, subset, weights, heterogeneity = c("iterative","unbiased"),TOL = 1E-6, echo = FALSE )
{
    cl <- match.call()  

    ## Get the appropriate column names or numbers.  
    form <- formula[[length(formula)]]
    require(nlme)
    val <- splitFormula(asOneSidedFormula(form[[3]]), sep = "+")    
    if("." %in% (years <- unlist(lapply(val, function(el) deparse(el[[2]])))))
        years <- grep("Y", names(x), perl = TRUE)

    ## Create a new data frame containing the appropriate contracts
    ## and years of experience.
    if (missing(subset))
        r <- TRUE
    else{
        e <- substitute(subset)
        r <- eval(e, data, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    ratios <- data[r, years]
      
    ## If weights are not specified, use equal weights as in
    ## Bühlmann's model.
    if (missing(weights))
    {
        if (any(is.na(ratios)))
            stop("missing values are not allowed in the matrix of ratios when the matrix of weights is not specified")
        weights <- array(1, dim(ratios))
        model <- "Bühlmann"
    }
    else
        #if (!identical(dim(weights), dim(ratios))) ## dimensions des weights? comment on donne ca?
            #stop
        model <- "Bühlmann-Straub"

    ## Check other bad arguments.
    if (ncol(ratios) < 2)
        stop("there must be at least one contract with at least two years of experience")
    if (nrow(ratios) < 2)
        stop("there must be more than one contract")
    if(!identical(which(is.na(ratios)), which(is.na(weights))))
        stop("missing values are not in the same positions in the matrix of weights and the matrix of ratios")

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

    if (identical(heterogeneity, "iterative"))
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

                    if (abs((at - at1)/at1) < TOL)
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

    ## Credibility premiums.
    P <- ratios.zw + cred * (ratios.w - ratios.zw)

    res <- list(model = model,
                premiums = P,
                individual = ratios.w,
                collective = ratios.zw,
                weights = weights.s,
                ncontracts = ncontracts,
                contracts = data[r,]$contract,
                s2 = s2,
                cred = cred,
                call = cl,
                unbiased = ac,
                iterative = at)
    class(res) <- "cm"
    res
}

print.cm <- function(x, ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Credibility premiums using the", x$model, "model: \n\n")
    sapply(paste("  Contract ",
                 format(x$contracts, justify = "right"),
                 ": ",
                 format(x$premiums),
                 sep = ""),
           cat, fill = TRUE)
    invisible(x)
}
