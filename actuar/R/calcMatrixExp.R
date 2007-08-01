### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the exponential matrix of exp( Q*u ) through
### diagonalisation of matrix Q and returns Pi_plus * exp( Q*u ) * 1_m
### where Pi_plus is the initial probability vector, Q the
### subintensity matrix and 1_m the one vector of R^m
###
### function used by ruinProb
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,


calcMatrixExp <- function(param)
{
    Q <- param$Q
    piplus <- param$piplus
    m <- length(piplus)
    ones <- rep(1,m)
    
    #compute the constant Ci through diagonalisation of Q
    resDiagQ <- eigen(Q)
    
    eigenValueQ <- resDiagQ$values
    P <- resDiagQ$vectors
    Pinv <- solve(P,diag(m))

    fctAux <- function(i) 
    {
        M<-array(0,c(m,m))
        M[i,i]<-1
        return(piplus %*% P %*% M %*% Pinv %*% ones)
    }
	
    Cste<-sapply(1:m,fctAux)

    psi <- function(u)
    {   
        sizeU <- length(u)
        t<-array(u,c(sizeU,1)) %*% array(eigenValueQ,c(1,m))
    
        # if there is a complex eigenvalue of Q, then there will be its conjugate in the spectrum of Q, so the terms
        # associated with this eigenvalue, and its conjugate in the ruin probability are real numbers
        c( Re( exp( t )  %*%  Cste ) )
    }
    
    return( psi )
}
