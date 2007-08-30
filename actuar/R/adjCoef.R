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

adjCoef <- function(mgf.claim, mgf.wait, premium, upper, h, reinsurance=NULL, from, to, n=101, toplot=FALSE)
{
    #check arguments
    if(missing(mgf.claim) && missing(h))
        stop("missing arguments: 'mgf.claim' or 'h' is needed")
    
    if(missing(upper))
        stop("missing argument 'upper', the upper bound of the m.g.f. of claim size distribution")

    if(!is.numeric(upper))
        stop("wrong argument 'upper', which must be a 'numeric'")

    if(!missing(h))
    {
        if(!is.function(h) || is.null(h))
            stop("wrong argument 'h', which must be a function")
    }
    
    if(!missing(mgf.claim))
    {
        if(!is.function(mgf.claim))
            stop("wrong argument 'mgf.claim', which must be a function")        
    
        if(!missing(mgf.wait))
        {
            if(!is.function(mgf.wait))
                stop("wrong argument 'mgf.wait', which must be a function")
        }
        else
            stop("missing argument 'mgf.wait'")
                
        if(missing(premium))
            stop("missing argument 'premium'")
    }  

    #reinsurance
    if(!is.null(reinsurance))
    {
        if(missing(from) || missing(to))
            stop("missing argument 'from' and/or 'to'")
        
        if(!is.numeric(from) || !is.numeric(to))
            stop("wrong argument 'from' and/or 'to', which must be a numeric")
        
        if(!is.function(premium) )
            stop("wrong argument 'premium', which must be a function when using reinsurance")

        if(!is.character(reinsurance) || ( reinsurance != "proportional" && reinsurance != "excess-of-loss") )
            stop("wrong argument 'reinsurance'")
    }

    #no reinsurance
    if(is.null(reinsurance) && !is.numeric(premium))
        stop("wrong argument 'premium', which must be a numeric")
    
    if(is.null(reinsurance))
    {      
        
        if(!missing(mgf.claim)) h <- function(r) mgf.claim(r)*mgf.wait(-r*premium)
 
        
    #compute the adjustment coefficient the unique positive root of h(r)=1
    #with h(r) = M_X(r) * M_W(-r*premRate) in the case of independence
    #where M_? stands for the moment generating function
        
        eqLundberg <-function(r)
        {
                                        #  cat("h ",h(r),"r ",r,"\n")
            ( h(r) - 1 )^2
        }

        interval <- c( 0, .9*upper )

    #minimisation through 'optimize'
        adjustcoeff <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))$minimum

        if(toplot)
        {
            if(eqLundberg(.3*upper) >1)
                x<-seq(0,.3*upper,length.out=100)
            else
                x<-seq(0,.6*upper,length.out=100)
            
            plot(x,sapply(x,h),type="l",col="blue",main="h",xlab="r",ylab="h(r)")
            lines(x,x*0+1,type="l")
        }

        comment(adjustcoeff) <- "no reinsurance"
        class(adjustcoeff) <- c("adjCoef", class(adjustcoeff))
        
        return(adjustcoeff)
    }

  
    #reinsurance
    if(!is.null(reinsurance))
    {
        retentionVect <- seq(from, to, length.out = n)

        #use to know if the Lundberg is defined by 'h' or mgf.{claim or wait}
        typeH <- !missing(mgf.claim)
        
        # function ret -> R(ret) 
        adjustCoefFunction <- function(ret)
        {
           if(typeH)
           {
               h <- function(r) mgf.claim(r,ret)*mgf.wait(-r*premium(ret))
               eqLundberg <-function(r)
                {
                    ( h(r) - 1 )^2
                }
           }
           else
           {    
               eqLundberg <-function(r)
               {
                   ( h(r,ret) -1 )^2
               }
           }
        
            interval <- c( 0, .9*upper )

            #minimisation through 'optimize'
            res <- optimize( eqLundberg, interval, tol = min(.Machine$double.eps^0.25, 10^-6))

            return(res$minimum)
        }
        
        vectAdjustCoef <- sapply( retentionVect, adjustCoefFunction)

        fun <- approxfun(retentionVect, vectAdjustCoef, rule=2, method = "linear")

        comment(fun) <- reinsurance
        class(fun) <- c("adjCoef", class(fun))
        attr(fun, "x") <- retentionVect
        attr(fun, "y") <- vectAdjustCoef

        return( fun )
    }
    
}

plot.adjCoef <- function(x, type = "l", ...)
{
    if(comment(x) %in% "no reinsurance")
    {
        warning("'x' is just a point not a function")
    }
    else
    {
        if(comment(x) %in% "proportional")
        {
            main <- "Adjustment Coefficient [Proportional reinsurance]"
            xlab <- "a [retention rate]"
            ylab <- "R(a)"
        }
        if(comment(x) %in% "excess-of-loss")
        {
            main <- "Adjustment Coefficient [Excess of loss reinsurance]"
            xlab <- "L [retention limit]"
            ylab <- "R(L)"
        }
        
        plot(attr(x, "x"), attr(x, "y"), main = main, xlab = xlab, ylab = ylab, type = type, ...)
    }
}

lines.adjCoef <- function(x, type = "l", ...)
{
    if(comment(x) %in% "no reinsurance")
    {
        warning("'x' is just a point not a function")
    }
    else
        lines(attr(x, "x"), attr(x, "y"), type = type, ...)
    
}


print.adjCoef <- function(x, ...)
{
    if(comment(x) %in% "no reinsurance")
    {
        cat("\nAdjustment coefficient\n")
        cat("value :", x, "\n")
    }
    if(comment(x) %in% "proportional")
    {
        cat("\nAdjustment coefficient - Proportional reinsurance\n")
        cat("\ngrid of retention rate a :", attr(x, "x"), "\n")
        cat("values at these points R(a) :", attr(x, "y"), "\n")
    }
    if(comment(x) %in% "excess-of-loss")
    {
        cat("\nAdjustment coefficient - Excess of loss reinsurance\n")
        cat("\ngrid of retention limit L :", attr(x, "x"), "\n")
        cat("values at these points R(L) :", attr(x, "y"), "\n")
    }
}
