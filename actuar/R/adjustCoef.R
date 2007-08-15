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

adjustCoef <- function(mgfClaim=NULL, mgfWaitTime=NULL, premRate=NULL, upper, h=NULL, toplot=FALSE)
{
    if(is.null(mgfClaim) && is.null(h))
        stop("missing argument: 'adjustCoef' needs either 'mgfClaim' or 'h' arguments")

    if(!is.null(mgfClaim) && !is.null(h))
        stop("too many arguments: 'adjustCoef' needs either 'mgfClaim' or 'h' arguments")

    if(is.null(upper))
        stop("missing argument 'upper', the upper bound of the m.g.f. of claim size distribution")

    if(!is.function(mgfClaim) && !is.null(mgfClaim) )
        stop("wrong argument 'mgfClaim', which must be a function")

    if(!is.function(mgfWaitTime) && !is.null(mgfWaitTime))
        stop("wrong argument 'mgfWaitTime', which must be a function")
              
    if(!is.function(h) && !is.null(h))
        stop("wrong argument 'h', which must be a function")

    if(!is.numeric(upper))
        stop("wrong argument 'h', which must be a 'numeric'")

    if(is.null(premRate) && !is.null(mgfClaim))
        stop("missing argument 'premRate', the premium rate")
   
    
    if(!is.null(mgfClaim))
    {
        if(!is.null(mgfWaitTime))
            h <- function(r) mgfClaim(r)*mgfWaitTime(-r*premRate)
        else
            h <- function(r) mgfClaim(r)*ratePoissonProcess/(ratePoissonProcess+r*premRate)
    }
           

    #compute the adjustment coefficient the unique positive root of h(r)=1
    #with h(r) = M_X(r) * M_W(-r*premRate) in the case of independence
    #where M_ stands for the moment generating function

    eqLundberg <-function(r)
        {
            cat("h ",h(r),"r ",r,"\n")
            ( h(r) - 1 )^2
        }

    interval <- c( 0, .9*upper )

    #minimisation through 'optimize'
    res <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))

    print(res)

    print(eqLundberg(.3*upper))
    
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

