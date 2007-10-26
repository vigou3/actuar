### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the subintensity matrix and the initial probability vector
### of the phase-type distribution parameters in the Sparre Andersen
### model.
###
### function used by ruinProb
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,

ruinSparreAndersen <- function(probpi,T,m,premRate,probnu,S,n)
{
    # claim sizes ~ PH(pi,T,m)
    # inter-occurence times ~ PH(nu,S,n)

    onesM <- rep(1,m)
    onesN <- rep(1,n)
    idM <- diag(m)
    idN <- diag(n)

    tzero <- -T %*% onesM
    szero <- -S %*% onesN
    
    matNu <- probnu %*% idN
    matPi <- probpi %*% idM

    
    
    # the contractant* (?) function phi, such that S verifies the point fixed equation Q=phi(Q)
	# the contractancy ensures (1) point fixed algorithm will converge to the (unique) solution
	# and from any starting value, (2) the convergence speed is exponential.
	# TODO : search in Berman & Plemmons (1979)** if phi is constractant.
	#* if contractant, we have it exists 0<k<1 such that |phi(x)-phi(y)| < k*|x-y|. 
	#** A. Berman & J. Plemmons (1979), Nonnegative matrices in the mathematical sciences, 
	#Academic Press, NY
    phi<-function(K)
    {
        # K "+" S with the rule for the kronecker sum for square matrix
        KplusS <- K %x% idN + idM %x% S
            
        Ahat <- -(idM %x% matNu) %*% solve( KplusS ,diag(m*n) ) %*% (idM %x% szero)
        
        return(T+tzero %*% matPi %*% Ahat)
    }

    #compute Q
    Qnew <- T
    tol <- 1
    preventLooping <- 0
    
    while(tol > 10^(-9) && preventLooping<75)
    {
        preventLooping <- preventLooping+1

        Qold <- Qnew
        Qnew <- phi(Qold)

        #infinite norm for matrix
        tol <- max(rowSums(abs(Qnew-Qold)))        
    }
    
    piplus <- onesM %*% (Qnew-T) /(sum(tzero)*premRate)
    res <- list(Q=Qnew,piplus=piplus)
    
    return(res)
}
