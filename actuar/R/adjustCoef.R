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
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,

adjustCoef <- function(model="CramerLundberg",marginal="exp",param,premRate,toplot=FALSE)
{
    #the user specifies directly the 'h' function
    if(is.function(model))
    {
        RBound <- param$RBound
        h <- model
        findRoot <- TRUE
    }
    else
    {
    #creation of the 'h' function
        if(model == "CramerLundberg")
        {
            findRoot <- FALSE
        #exponential claim sizes E(beta)
            if(marginal == "exp")
            {
                beta <- param$beta
                lambda <- param$lambda
                return( beta-lambda/premRate )
            }
        #gamma claim sizes G(alpha,beta)
            else if(marginal == "gamma")
            {
                lambda <- param$lambda
                alpha <- param$alpha
                beta <- param$beta
                
                if(length(beta) != 1) return()
                
                h<-function(r) (beta/(beta-r))^alpha*lambda/(lambda+r*premRate)
                RBound <- beta
                findRoot <- TRUE
            }
        #claim sizes are a mixture of exponential MixExp(weight,beta)
            else if(marginal == "mixexp")
            {
                lambda <- param$lambda
                weight <- param$weight
                beta <- param$beta
                
                if(length(beta) != length(weight)) return()
                
                h<-function(r)
                    sum(weight*beta/(beta-r)) * lambda/(lambda+r*premRate)
                
                RBound <- min(beta)
                findRoot <- TRUE            
            }
        #claim sizes are a mixture of gamma MixGam(weight,alpha,beta)
            else if(marginal == "mixgamma")
            {
                lambda <- param$lambda
                weight <- param$weight
                beta <- param$beta
                alpha <- param$alpha
                
                if(length(beta) != length(weight) || length(beta) != length(alpha)) return()

                h<-function(r)
                    sum( weight*(beta/(beta-r))^alpha ) * lambda/(lambda+r*premRate)
                
                RBound <- min(beta)
                findRoot <- TRUE
            }
        #otherwise phase type claim size PH(pi,T,m)
            else if(marginal == "phasetype")
            {
                lambda <- param$lambda
                probpi <- param$pi
                matT <- param$T
                m <- param$m

                ones <- rep(1,m)       
                tzero <- -matT %*% ones
                idM <- diag(m)
                                
                mgfClaimSize<-function(s)                          
                    probpi %*% solve(-s*idM - matT,idM) %*% tzero
                

                h<-function(r)
                    mgfClaimSize(r) * lambda/(lambda+r*premRate)
                
                RBound <- min(-diag(matT))
                findRoot <- TRUE
                
            }
        }
    #creation of the 'h' function
        else if(model == "SparreAndersen")
        {
        #exponential claim sizes E(beta)
        #gamma inter-occurence times G(eta,lambda)
            if(marginal == "expgamma")
            {
                beta <- param$beta
                lambda <- param$lambda
                eta <- param$eta

                h <- function(r)
                    beta/(beta-r)*(lambda/(lambda+r*premRate))^eta
                
                findRoot <- TRUE
                RBound <- beta
            }
        #gamma claim sizes G(eta,beta)
        #gamma inter-occurence times G(eta,lambda)
            else if(marginal == "gammagamma")
            {
                beta <- param$beta
                alpha <- param$alpha
                lambda <- param$lambda
                eta <- param$eta
                            
                h <- function(r)
                    (beta/(beta-r))^alpha*(lambda/(lambda+r*premRate))^eta
            
                findRoot <- TRUE
                RBound <- beta
            }
        #phase type claim sizes PH(pi,T,m)
        #phase type inter-occurence times PH(nu,S,n)
            else if (marginal == "phasetype")
            {
                nu <- param$nu
                matS <- param$S
                n <- param$n
                
                probpi <- param$pi
                matT <- param$T
                m <- param$m
                
                onesM <- rep(1,m)       
                tzero <- -matT %*% onesM
                idM <- diag(m)
                
                onesN <- rep(1,n)       
                szero <- -matS %*% onesN
                idN <- diag(n)
                
                mgfClaimSize<-function(s)             
                    probpi %*% solve(-s*idM - matT,idM) %*% tzero
                
                mgfInterOccurTime<-function(s)           
                    nu %*% solve(-s*idN - matS,idN) %*% szero
                
                h<-function(r)
                    mgfClaimSize(r) *mgfInterOccurTime(-r*premRate)
                                
                RBound <- min(-diag(matT))
                findRoot <- TRUE
            }

        
        }
        else return()
    }
    
    

    #compute the adjustment coefficient the unique positive root of h(r)=1
    #with h(r) = M_X(r) * M_W(-r*premRate) in the case of independence
    #where M_ stands for the moment generating function
    if(findRoot)
    {
        eqLundberg <-function(r)
            (h(r) - 1)^2

        interval <- c( 0, .9*RBound )
        
        res <- optimize( eqLundberg ,interval)
        

        if(toplot)
        {
            if(eqLundberg(.3*RBound) >1)
                x<-seq(0,.3*RBound,length.out=100)
            else
                x<-seq(0,.6*RBound,length.out=100)
            
            plot(x,sapply(x,h),type="l",col="blue",main="h",xlab="r",ylab="h(r)")
            lines(x,x*0+1,type="l")
        }
        
        return(res$minimum)
    }
    else
        return()
}

