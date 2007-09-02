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

adjCoef <- function(mgf.claim, mgf.wait, premium, upper, h, reinsurance=c("none","proportional","excess-of-loss"), from, to, n=101, toplot=FALSE)
{
    reinsurance <- match.arg(reinsurance)
    
    ## check arguments
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
    if(reinsurance != "none")
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
    if(reinsurance == "none" && !is.numeric(premium))
        stop("wrong argument 'premium', which must be a numeric")

    
    ## compute the adjustment coefficient

    #non reinsurance
    if(reinsurance == "none")
    {      
        
        if(!missing(mgf.claim)) h <- function(r) mgf.claim(r)*mgf.wait(-r*premium)
 
        
    #compute the adjustment coefficient the unique positive root of h(r)=1
    #with h(r) = M_X(r) * M_W(-r*premRate) in the case of independence
    #where M_? stands for the moment generating function
        
        eqLundberg <-function(r) ( h(r) - 1 )^2

        interval <- c( 0, upper-.Machine$double.eps )

    #minimisation through 'optimize'
    # le 'return' est obligatoire car sinon R continue à exéctuer le code
    # du coup, il renvoie NULL
        return(optimize( eqLundberg, interval, tol = sqrt(.Machine$double.eps))$minimum)
    }

  
    #reinsurance
    if(reinsurance != "none")
    {
        retentionVect <- seq(from, to, length.out = n)
        
        #if the mgf.{claim and wait} is defined
        if(!missing(mgf.claim))
            h <- function(r,ret) mgf.claim(r,ret)*mgf.wait(-r*premium(ret))
        #otherwise 'h' is defined by the user
        
        # function ret -> R(ret) 
        adjustCoefFunction <- function(ret)
        {
            eqLundberg <-function(r) ( h(r,ret) -1 )^2           
        
            interval <- c( 0, upper-.Machine$double.eps )

            #minimisation through 'optimize'
            optimize( eqLundberg, interval, tol = sqrt(.Machine$double.eps))$minimum
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

plot.adjCoef <- function(x, type = "l", add = FALSE, ...)
{
    if("function" %in% class(x))
    {
        if(comment(x) == "proportional")
        {
            main <- "Adjustment Coefficient [Proportional reinsurance]"
            xlab <- "a [retention rate]"
            ylab <- "R(a)"
        }
        else
        {
            if(comment(x) != "excess-of-loss")
                stop("internal error in adjCoef, wrong comment for an 'adjCoef' object")
            
            main <- "Adjustment Coefficient [Excess of loss reinsurance]"
            xlab <- "L [retention limit]"
            ylab <- "R(L)"
        }

        #Warning messages:
        #1: "add" n'est pas un paramètre graphique in: plot.window(xlim, ylim, log, asp, ...) 
        #2: "add" n'est pas un paramètre graphique in: plot.xy(xy, type, pch, lty, col, bg, cex, lwd, ...) 
        #3: "add" n'est pas un paramètre graphique in: axis(side, at, as.graphicsAnnot(labels), tick, line, pos, outer,  
        #4: "add" n'est pas un paramètre graphique in: axis(side, at, as.graphicsAnnot(labels), tick, line, pos, outer,  
        #5: "add" n'est pas un paramètre graphique in: box(which = which, lty = lty, ...) 
        #6: "add" n'est pas un paramètre graphique in: title(main, sub, xlab, ylab, line, outer, ...)
                                                   
        if(add) 
            lines(attr(x, "x"), attr(x, "y"), main=main, xlab=xlab, ylab=ylab, type=type,...)
        else
            plot(attr(x, "x"), attr(x, "y"), main=main, xlab=xlab, ylab=ylab, type=type,...)
    }
}

lines.adjCoef <- function(x, type = "l", ...)
{
    if("function" %in% class(x))
        lines(attr(x, "x"), attr(x, "y"), type=type, add=TRUE, ...)    
}

