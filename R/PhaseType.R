### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {d,p,r,m,mgf}ph functions to compute
### characteristics of Phase-type distributions with cumulative
### distribution function
###
###       Pr[X <= x] = 1 - pi %*% exp(Tx) %*% 1m,
###
### where 'pi' is the initial probability vector, 'T' is the
### subintensity matrix and 1m is 1-vector of R^m.
###
### See Bladt (2005).
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dphtype <- function (x, prob, transition, log = FALSE)
    .External("do_dpqph", "dph", x, prob, transition, log)

pphtype <- function (x, prob, transition, lower.tail = TRUE, log.p = FALSE)
    .External("do_dpqph", "pph", x, prob, transition, lower.tail, log.p)

rphtype <- function (n, prob, transition, log.p = FALSE)
    .External("do_randomph", "rph", x, prob, transition)

mphtype <- function (x, prob, transition)
    .External("do_dpqph", "mph", x, prob, transition, FALSE)

mgfphtype <- function (x, prob, transition, log = FALSE)
    .External("do_dpqph", "mgfph", x, prob, transition, log)


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

    survivalProb <- .Call("calcMatExpGen", q, pi, T, onesM)

    if(!log.p)
    {
        if(lower.tail)
            return(1 - survivalProb)
        else
            return(survivalProb)
    }
    else
    {
        if(lower.tail)
            return( log(1 - survivalProb) )
        else
            return( log(survivalProb) )
    }
}

qphasetype <- function(p, pi, T, m, lower.tail = TRUE, log.p = FALSE)
    return (NULL)
#???

rphasetype <- function(n, pi, T, m)
{
    if(length( dim(T) ) != 2)
    {
        stop("T is not a 'mxm' matrix")
    }
    if(dim(T)[1] != dim(T)[2])
    {
        stop("T is not a 'mxm' matrix")
    }

    tzero <- c(-T %*% rep(1,m))
    .Call("randphasetype", n, pi, T, tzero )
}


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

    routine <- function(x)
    {
        if(x == 0)
            return(1)
        if(x >= 1)
        {
        #integer order
            if(as.integer(x) == x)
            {
                TpowN <- T

                factN <- 1
                if(x > 1)
                { #prod(T==diag(diag(T)))
                    for(k in 2:x)
                    {
                        factN <- factN * k
                        TpowN <- TpowN %*% T
                        #cat("-",k)
                        #print(TpowN)
                    }
                #cat("\nn!",factN,"\n")
                }
                return( factN * (-1)^x * pi %*% solve(TpowN, diag(m)) %*% onesM )
            }
            else
            {
                stop("non integer orders are not supported")
            }
        }
    }
    return(sapply(order, routine))
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

     routine <- function(x)
     {
         options(show.error.messages = FALSE)
         #test the matrix inversion
         #if failed, testDiag is an invisible object of class 'try-error'
         testDiag <- try( solve(-x*diag(m) - T , diag(m) ) )
         options(show.error.messages = TRUE) #revert to default

         if(class(testDiag) == "try-error")
             return( NaN )
         else
             return( pi %*% solve(-x*diag(m) - T , diag(m) ) %*% tzero )
     }


     if(!log)
         return( sapply(x, routine) )
     else
         return( log( sapply(x, routine) ) )
}
