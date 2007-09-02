### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the exponential matrix of exp( Q*u ) through
### diagonalisation of matrix Q and returns Pi_plus * exp( Q*u ) * 1_m
### where Pi_plus is the initial probability vector, Q the
### subintensity matrix and 1_m the one vector of R^m
###
### function used by ruin
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,


calcMatrixExp <- function(param)
{
    Q <- param$Q
    piplus <- param$piplus
    m <- length(piplus)
    ones <- rep(1,m)

    # try diagonalisation of Q
    options(show.error.messages = FALSE)
    #test the diagonalisation of Q
    #if failed, testDiag is an invisible object of class 'try-error'
    testDiag <- try( solve( eigen(Q)$vectors , diag(m) ) )
    options(show.error.messages = TRUE) #revert to default
            
    if(class(testDiag) == "try-error")
    {
        ## compute matrix exponential in the general case

        #old error
                                        #stop("\n\t*** actuar internal error : Q is non diagonalisable ***\n\t",geterrmessage()) #stop execution of ruinProb
        psi <- function(u)
            .Call("calcMatExpGen", piplus, u, Q, ones)

        return( psi )
    }
    else
    {
        ## compute matrix exponential through diagonalisation of Q
        return( calcMatrixExpDiag(param) )
    }
}

