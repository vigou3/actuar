### ===== actuar: an R package for Actuarial Science =====
###
### 'Predict' method for Generalized Hierarchical Credibility Models
###
### Predicted values of claim ratios (credibility premiums)
### based on hierarchical credibility model object.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot

predict.ghcm <- function(object, ...)
{
    nLevels <- length(object$means)
    premiums <- vector("list", nLevels)
    cred <- rev(object$cred)
    premiums[[1]] <- as.vector(object$means[[nLevels]])
    aff <- rev(object$aff)
    for (i in 2:nLevels)
    {
        M <- rep(premiums[[i-1]], table(aff[[i-1]]))[order(aff[[i-1]])]
        cred <- as.vector(rev(object$cred)[[i-1]])
        means <- as.vector(rev(object$means)[[i]])        
        premiums[[i]] <- cred * means + (1 - cred) * M
        names(premiums[[i]]) <- unique(object$data[rev(object$levs)][[i]])
    }
    names(premiums) <- rev(object$levs)
    premiums
}
