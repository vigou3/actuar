### ===== actuar: an R package for Actuarial Science =====
###
### Ruin Theory
###
### Compute the adjustment coefficient R, the (strictly) positive
### root of the Lundberg equation :
###         h(r) = E[ e^(r X - r c W) ] = 1
### where X is the generic claim size random variable,
### W the inter-occurence time, and c the premium rate.
### If c does not respect the postive safety loading constraint
### E[X-cW]<0, R does not exist, the function returns 0,
### otherwise, it returns R>0.
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>
###

adjustCoef <- function(mgfClaim, mgfWaitTime, premRate, upper, h, retRate=1, retLimit=Inf, toplot=FALSE,...)
{
    #check arguments
    if(missing(mgfClaim) && missing(h))
        stop("missing arguments: 'adjustCoef' needs either 'mgfClaim' or 'h' arguments")
    
    if(!missing(mgfClaim) && !missing(h))
        stop("too many arguments: 'adjustCoef' needs either 'mgfClaim' or 'h' arguments")

    if(missing(upper))
        stop("missing argument 'upper', the upper bound of the m.g.f. of claim size distribution")

    if(!is.numeric(upper))
        stop("wrong argument 'upper', which must be a 'numeric'")

    if(!missing(h))
    {
        if(!is.function(h) || is.null(h))
            stop("wrong argument 'h', which must be a function")
    }
    
    if(!missing(mgfClaim))
    {
        if(!is.function(mgfClaim) || is.null(mgfClaim))
            stop("wrong argument 'h', which must be a function")        
    
        if(!missing(mgfWaitTime))
        {
            if(!is.function(mgfWaitTime) || is.null(mgfWaitTime))
                stop("wrong argument 'mgfWaitTime', which must be a function")
        }
        
        if(missing(mgfWaitTime) && missing(ratePoissonProcess))
            stop("missing argument 'ratePoissonProcess', the rate of the Poisson process")
        
        if(missing(premRate))
            stop("missing argument 'premRate', the premium rate")
    }  

    if(length(retRate) > 1 || length(retLimit) > 1)
    {
        #proportional reinsurance
        if(length(retRate) > 1 && !is.vector(retRate) )
            stop("wrong argument 'retRate', which is must be a 'vector'")
        
        if(length(retRate) > 1 && !is.function(premRate) )
            stop("wrong argument 'premRate', which must be a function when using proportional reinsurance")

        #excess of loss reinsurance
        if(length(retLimit) > 1 && !is.vector(retLimit) )
            stop("wrong argument 'retLimit', which is must be a 'vector'")
        
        if(length(retLimit) > 1 && !is.function(premRate) )
            stop("wrong argument 'premRate', which must be a function when using excess of loss reinsurance")
    }

    #no reinsurance
    if(length(retLimit) == 1 && length(retRate) == 1 && !is.numeric(premRate) && !missing(mgfClaim))
        stop("wrong argument 'premRate', which must be a numeric")
    
    if(retRate[1] == 1 && retLimit[1] == Inf)
    {      
        
        if(!missing(mgfClaim))
        {
            if(!missing(mgfWaitTime))
                h <- function(r) mgfClaim(r)*mgfWaitTime(-r*premRate)
            else
                h <- function(r) mgfClaim(r)*ratePoissonProcess/(ratePoissonProcess+r*premRate)
        }
        
        
    #compute the adjustment coefficient the unique positive root of h(r)=1
    #with h(r) = M_X(r) * M_W(-r*premRate) in the case of independence
    #where M_ stands for the moment generating function
        
        eqLundberg <-function(r)
        {
                                        #  cat("h ",h(r),"r ",r,"\n")
            ( h(r) - 1 )^2
        }

        interval <- c( 0, .9*upper )

    #minimisation through 'optimize'
        res <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))

        if(toplot)
        {
            if(eqLundberg(.3*upper) >1)
                x<-seq(0,.3*upper,length.out=100)
            else
                x<-seq(0,.6*upper,length.out=100)
            
            plot(x,sapply(x,h),type="l",col="blue",main="h",xlab="r",ylab="h(r)")
            lines(x,x*0+1,type="l")
        }
        
        return(res$minimum)
    }

    if(retRate[1] != 1 && retLimit[1] != Inf)
        stop("a mix of excess of loss and proportional reinsurance is not supported")

    #proportional reinsurance
    if(retRate[1] != 1)
    {
        typeH <- 0
        if(!missing(mgfClaim) && !missing(mgfWaitTime)) typeH <- 1
        if(!missing(mgfClaim) && missing(mgfWaitTime)) typeH <- 2
        if(!missing(h)) typeH <- 3

        if(typeH == 0) stop("internal error in adjustCoef")
        
        # function a -> R(a) 
        adjustCoefFunction <- function(a)
        {
           if(typeH == 1)
               h <- function(r) mgfClaim(r,a)*mgfWaitTime(-r*premRate(a))
           if(typeH == 2)
               h <- function(r) mgfClaim(r,a)*ratePoissonProcess/(ratePoissonProcess+r*premRate(a))                               
            

            if(typeH == 1 || typeH == 2)
            {
                eqLundberg <-function(r)
                {
                                        #  cat("h ",h(r),"r ",r,"\n")
                    ( h(r) - 1 )^2
                }
            }
            if(typeH == 3)
            {
                eqLundberg <-function(r)
                {
                    ( h(r,a) -1 )^2
                }
            }
            
            interval <- c( 0, .9*upper )

            #minimisation through 'optimize'
            res <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))

            return(res$minimum)
        }
        
        vectAdjustCoef <- sapply(retRate, adjustCoefFunction)

        return( approxfun(retRate, vectAdjustCoef, rule=2, method = "linear") )
    }
    
    
    #excess of loss reinsurance
    if(retLimit[1] != Inf)
    {
        typeH <- 0
        if(!missing(mgfClaim) && !missing(mgfWaitTime)) typeH <- 1
        if(!missing(mgfClaim) && missing(mgfWaitTime)) typeH <- 2
        if(!missing(h)) typeH <- 3

        if(typeH == 0) stop("internal error in adjustCoef")

        print(typeH)
        
        # function L -> R(L) 
        adjustCoefFunction <- function(L)
        {
           if(typeH == 1)
               h <- function(r) mgfClaim(r,L)*mgfWaitTime(-r*premRate(L))
           if(typeH == 2)
               h <- function(r) mgfClaim(r,L)*ratePoissonProcess/(ratePoissonProcess+r*premRate(L))                               
           

            if(typeH == 1 || typeH == 2)
            {
                eqLundberg <-function(r)
                {
                                        #  cat("h ",h(r),"r ",r,"\n")
                    ( h(r) - 1 )^2
                }
            }
            if(typeH == 3)
            {
                eqLundberg <-function(r)
                {
                    #cat("h ",h(r,L),"r ",r,"L ",L,"\n")
                    ( h(r,L) -1 )^2
                }
            }
            
            interval <- c( 0, .9*upper )

            #minimisation through 'optimize'
            res <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))
           
           if(toplot)
           {
               if(eqLundberg(.3*upper) >1)
                   x<-seq(0,.3*upper,length.out=100)
               else
                   x<-seq(0,.6*upper,length.out=100)
               
               plot(x,sapply(x,h,L),type="l",col="blue",main="h",xlab="r",ylab="h(r)")
               lines(x,x*0+1,type="l")
           }

            return(res$minimum)
        }
        print("blii")
        
        vectAdjustCoef <- sapply(retLimit, adjustCoefFunction)

        
        

        return( approxfun(retLimit, vectAdjustCoef, rule=2, method = "linear") )
    }
    
        
    
}

