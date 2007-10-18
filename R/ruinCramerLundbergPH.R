### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the subintensity matrix and the initial probability vector
### of the phase-type distribution parameters in the Cramér-Lundberg
### model.
###
### function used by ruinProb
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,


ruinCramerLundbergPH <- function(probpi,T,m,premRate,lambda)
{
        #check the sizes of the argument, and remove the case where claim size is exponential
    if(length(attributes(T)$dim) != 2 )
      return()

        #compute the matrix Q, sub-intensity matrix of the ruin probability
    Tinv <- solve(T,diag(m))
    ones <- rep(1,m)
	
    tzero <- -T %*% ones
    piplus <- -lambda/premRate*probpi %*% Tinv
    
    Q <- T+tzero%*%piplus
            
    return(list(Q=Q,piplus=piplus))   	
}
