### ===== actuar: an R package for Actuarial Science =====
###
### Fitting Credibility Models:
###
### Fit a credibility model in the formulation of variance components 
### as described in Dannenburg, Kaas and Goovaerts (1996). Models supported
### are part of a generalized hierarchical credibility theory as introduced
### in Dannenburg (1995).
### 
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot

cm <- function(formula, data = NULL, subset, weights, TOL=1E-6, echo=FALSE)
{
    cl <- match.call()
    ## Create a new data frame containing the appropriate contracts
    ## and years of experience.
    
    ## Extract the name of each variable.
    vars <- all.vars(formula)
    years <- all.vars(asOneSidedFormula(formula[[2]][[3]]))
    levs <- all.vars(lhs <- asOneSidedFormula(formula[[2]][[2]]))
    pfstruct <- c(rev(rownames(attr(terms(lhs),"factors"))), "pf")

    if (missing(data))
    {
        ## Stop if vectors of data in 'formula' are of different lengths.
        ## There must be a nicer way...
        l <- sapply(vars, function(x) length(get(x)))
        ll <- sapply(l, function(x) x == l[1])
        if (!all(ll))
            stop("all objects in 'formula' must be of equal length")
        if ("." %in% years)
            stop("'.' operator is meaningless if 'data' is not specified")
        data <- data.frame(sapply(vars, get))        
    }
    else
    {
        if (!inherits(data, "data.frame"))
            stop("'data' must be a 'data frame'")
        if (!all(vars %in% c(names(data), ".")))
            stop("variables in 'formula' must be names from 'data' or existing 'R' objects") 
        if ("." %in% years)
            years <- suppressWarnings(names(data)[names(data) != pfstruct])       
    }
    ## Bind a column of '1's representing the affiliation to the one portfolio.
    ## Used to symmetrize the further calculations.    
    data <- cbind(pf  = rep(1, nrow(data)), data[c(levs, years)])

    ## If weights are not specified, use equal weights.
    if (missing(weights)) 
    {
        if (any(is.na(data[years])))
            stop("missing values of ratios are not allowed when the matrix of weights is not specified")
        weights <- data
        weights[names(weights) != pfstruct] <- sapply(weights[years], function(x) rep(1, length(x)))
    }
    else
    {
        if (ncol(data[years]) < 2)
            stop("there must be at least one contract with at least two years of experience")
        if (nrow(data) < 2)
            stop("there must be more than one contract")        
        weights <- cbind(data[pfstruct], weights)
        names(weights)[names(weights) != pfstruct] <- years
        if (!identical(which(is.na(data)), which(is.na(weights))))
            stop("missing values are inconsistent in 'ratios'/'weights' data")            
    }
    
    nlevels <- length(pfstruct) - 1    
    nstruct <- c(length(years), sapply(pfstruct, function(x) length(unique(data[[x]]))))

    ## s^2
    ind.weight <- rowSums(weights[years], na.rm = TRUE)
    ind.means <- rowSums(data[years]*weights[years], na.rm = TRUE)/ind.weight
    s2num <- sum(weights[years]*(data[years]-ind.means)^2)
    denoms <- numeric(nlevels + 1)
    
    for (i in 2:(length(denoms))) denoms[i] <- nstruct[i] - nstruct[i+1]
    s2 <- s2num/(denoms[1] <- length(na.omit(unlist(data[years]))) - nstruct[[2]])

    ## Create vectors for values to be outputted.
    
    param <- rep(s2, nlevels + 1)        # The structure parameters
    cred <- vector("list", nlevels)      # The credibility factors
    w. <- vector("list", nlevels)        # The credibility weights
    M <- vector("list", nlevels)         # The individual and collective estimators  
     
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
        for (i in 2:(nlevels + 1))
        {
            cred[[i-1]] <- 1/(1 + param[i-1]/(param[i] * weight))
            aff <- tapply(data[[pfstruct[i]]], data[[pfstruct[i-1]]], function(x) unique(x)) 
            weight. <- tapply(cred[[i-1]], aff, sum)
            means. <- tapply(cred[[i-1]] * means, aff, sum) / weight.
            param[i] <- sum(cred[[i-1]] *
                            (means - rep(means., table(aff))[order(aff)])^2) / denoms[i]
            
            w.[[i-1]] <- weight
            weight <- weight.

            M[[i-1]] <- means
            means <- means.
        }
        p <- ifelse(any(param <= TOL), which(param > TOL), TRUE)
        if (max(abs((param[p] - paramt[p])/paramt[p])) < TOL) 
                break        
    }
    names(param) <- c("s2", letters[1:nlevels])
    list(param = param, weights = w., means = M, cred = cred)
}
                                        
    
if (missing(subset))
    r <- TRUE
else{
    e <- substitute(subset)
    r <- eval(e, data, parent.frame())
    if (!is.logical(r)) 
        stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
    if (length(unique(data[[struct]])) < length(r[r]))    ####### probleme... a voir
        stop("Hierarchical conflict in 'formula'/'subset'")        
}
if (length(unique(data[[struct]])) != length(data[[struct]]))
    data <- t(as.data.frame(lapply(by(data, data[struct],
                                      match.fun(quote(subset)),
                                      select = years),
                                   colSums)))
else 







        
            
                              

        
