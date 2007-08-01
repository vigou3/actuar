### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the infinite time ruin probability in the model
### of Cramér-Lundberg or Sparre Andersen, using the results
### of Gerber-Dufresnes (1988) and Asmussen-Rolski (1991).
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,

ruinProb <-function(model="CramerLundberg",param,premRate)
{
    
    ## Check arguments (TODO)
    
    if(model == "CramerLundberg")
    {
        lambda <- param$lambda
        probpi <- param$pi
        T <- param$T
        m <- param$m
        
        ## exponential claim size
        if(m == 1)
        {            
            beta <- T
            return( function(u) lambda / (beta*premRate) *exp(-(beta-lambda/premRate)*u) )
        }

        gotoPhaseType <- FALSE
        
        ## mixture of 2 exponentials - 'Gerber Dufresnes method'
        if(m == 2)
        {
            if(T[2,1] == 0 && T[1,2] == 0)
            {
                beta <- -diag(T)
              
                ## the discriminant of the second order equation
                delta <- 1+(premRate/lambda)^2*(beta[1]-beta[2])^2 + 2*premRate/lambda*(sum(beta)-2*sum(beta*probpi) )

                if(delta >0)
                {
                    r1 <- (-lambda+premRate*(sum(beta))-lambda*sqrt(delta))/(2*premRate)
                    r2 <- (-lambda+premRate*(sum(beta))+lambda*sqrt(delta))/(2*premRate)
                    
                    C1 <- r2/(r2-r1)*(beta[1]-r1)/beta[1]*(beta[2]-r1)/beta[2]
                    C2 <- r1/(r1-r2)*(beta[2]-r2)/beta[2]*(beta[1]-r2)/beta[1]
                   
                    return( function(u) C1*exp(-r1*u)+ C2*exp(-r2*u) )
                }
                else gotoPhaseType <- TRUE

            }
            gotoPhaseType <- TRUE
        }

        ## claim sizes are phase type - 'Asmussen method'
        if(m >= 3 || gotoPhaseType)
        {
            
            ## compute the sub intensity matrix and the initial distribution
            ## of the phase type distribution of the ruin probability
            resQPiplus <- ruinCramerLundbergPH(probpi,T,m,premRate,lambda)

            ## check the difference of psi(1) with package matrix and diagonalisation of Q
            Q <- resQPiplus$Q
            piplus <- resQPiplus$piplus

            resDiagPsi <- calcMatrixExp(resQPiplus)
            resDiag <- resDiagPsi(1)
            
            resPackageMatrix <- piplus %*% as.matrix(expm(Matrix(Q))) %*% rep(1,length(piplus))
            
            error <- max(rowSums(abs(resDiag-resPackageMatrix)))

            
            if(error < 10^-6)
            {
                ## compute matrix exponential through diagonalisation of Q                
                return(calcMatrixExp(resQPiplus))
            }
            else            
                return( function(u) piplus %*% as.matrix(expm(Matrix(Q)*u)) %*% rep(1,length(piplus)) )                        
        }
    }

    if(model == "SparreAndersen")
    {
        probpi <- param$pi
        T <- param$T
        m <- param$m

        probnu <- param$nu
        S <- param$S
        n <- param$n
        
        ## compute the sub intensity matrix and the initial distribution
        ## of the phase type distribution of the ruin probability
        resQPiplus <- ruinSparreAndersen(probpi,T,m,premRate,probnu,S,n)

        ## check the difference of psi(1) with package matrix and diagonalisation of Q
        Q <- resQPiplus$Q
        piplus <- resQPiplus$piplus

        resDiagPsi <- calcMatrixExp(resQPiplus)
        resDiag <- resDiagPsi(1)
            
        resPackageMatrix <- piplus %*% as.matrix(expm(Matrix(Q))) %*% rep(1,length(piplus))
        
        error <- max(rowSums(abs(resDiag-resPackageMatrix)))
        
        if(error < 10^-6)
        {
            ## compute matrix exponential through diagonalisation of Q            
            return(calcMatrixExp(resQPiplus))
        }
        else
            return( function(u) piplus %*% as.matrix(expm(Matrix(Q)*u)) %*% rep(1,length(piplus)) )                                
    }
}
