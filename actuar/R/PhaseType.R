### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,mgf}phasetype functions to compute
### characteristics of the Phase-type distribution. Its distribution function
### is 
###       P(X <= x) = 1-pi %*% exp(T) %*% 1m
### where pi is the initial probability vector, T the subintensity matrix and
### 1m is 1-vector of R^m
###
### See Bladt(2005).
### 
###
### AUTHORS:  Mathieu Pigeon, Christophe Dutang and
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

dphasetype <- function (x, pi, T, m, log = FALSE)
{
    if(length( dim(T) ) != 2)
    {
        stop("T is not a 'mxm' matrix")
    }
    if(dim(T)[1] != dim(T)[2])
    {
        stop("T is not a 'mxm' matrix")
    }    
    
    onesM <- rep(1,m)
    tzero <- c(-T %*% onesM)

    if(!log)
        return( .Call("calcMatExpGen", x, pi, T, tzero) )
    else
        return( log( .Call("calcMatExpGen", x, pi, T, tzero) ) )
}

pphasetype <- function(q, pi, T, m, lower.tail = TRUE, log.p = FALSE)
{
    if(length( dim(T) ) != 2)
    {
        stop("T is not a 'mxm' matrix")
    }
    if(dim(T)[1] != dim(T)[2])
    {
        stop("T is not a 'mxm' matrix")
    }    
    
    onesM <- rep(1,m)

    survivalProb <- .Call("calcMatExpGen", x, pi, T, onesM)

    if(!log)
        return(1 - survivalProb)
    else
        return( log(1 - survivalProb) )
}

qphasetype <- function(p, pi, T, m, lower.tail = TRUE, log.p = FALSE)
    return (NULL)
#???
    
rphasetype <- function(n, pi, T, m)
    .Call("randphasetype", n, pi, T, -T %*% rep(1,m) )


mphasetype <- function(order, pi, T, m)
{
    if(length( dim(T) ) != 2)
    {
        stop("T is not a 'mxm' matrix")
    }
    if(dim(T)[1] != dim(T)[2])
    {
        stop("T is not a 'mxm' matrix")
    }    
    
    onesM <- rep(1,m)
    
    if(order == 0)
        return(1)
    if(order >= 1)
    {
        #integer order
        if(as.integer(order) == order)
        {
            TpowN <- T
            n <- order-1
            factN <- 1
            for(k in 1:n)
            {
                factN <- factN * k
                TpowN <- TpowN %*% TpowN
            }
            return( factN * (-1)^order * pi %*% solve(TpowN, diag(m)) %*% onesM )
        }
        else
        {
            stop("non integer order is not supported")
        }
    }
    
}

mgfphasetype <- function(x, pi, T, m, log = FALSE)
{
     if(length( dim(T) ) != 2)
    {
        stop("T is not a 'mxm' matrix")
    }
    if(dim(T)[1] != dim(T)[2])
    {
        stop("T is not a 'mxm' matrix")
    }    
    
    onesM <- rep(1,m)
    tzero <- c(-T %*% onesM)
     
    temp <- solve(-x*diag(m) - T , diag(m) )

     if(!log)
         return( pi %*% temp %*% tzero )
     else
         return( log( pi %*% temp %*% tzero ) )
}
    

