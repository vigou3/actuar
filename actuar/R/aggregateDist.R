aggregateDist <- function(method = c("normal", "np2", "simulation", "recursive", "exact"),
                          model.sev, model.freq, moments = c(mean = 0, var= 1, skewness = NULL),
                          x.scale = 1, n, p0, TOL = 1e-06, echo = FALSE, ...)
{
    
    ## The method used essentially tells which function should be
    ## called for the calculation of the aggregate claims
    ## distribution.
  
    method <- match.arg(method)
    
    
    if (method == "normal"){
       
        ## An error message is issued if the number of moments listed
        ## is not appropriate regarding the method. However it is the
        ## user's responsability to list the moments in the correct
        ## order since the vector is not required to be named.

             
        if (length(moments) != 2) stop("'normal' method requires the first TWO moments of the distribution")
        return(normal(moments[1], moments[2]))}
        
    if (method == "np2"){
        if (length(moments) != 3) stop("'np2' method requires the first THREE moments of the distribution")
        return(np2(moments[1], moments[2], moments[3]))}

    
    if (method == "simulation") return(simS(n, model.freq, model.sev))

    ## If 'model.sev' or 'model.freq' are vectors of probabilities,
    ## they are directly passed on to the subfunction. If they are
    ## expressed as parameterized distributions, they are discretized
    ## before being passed on.

    if (class(model.sev) == "numeric") fx <- model.sev
       
    else
    {
        psev <- match.fun(paste("p", model.sev$dist, sep = ""))
        qsev <- match.fun(paste("q", model.sev$dist, sep = ""))
        qsevpar <- c(p = 1 - 10e-04*TOL, model.sev$par)
        formals(qsev)[names(qsevpar)]  <- qsevpar
        formals(psev)[names(model.sev$par)] <- model.sev$par
        Fx <- psev(seq(0, qsev()))  
        fx <- c(0, diff(Fx))
    }
    if (method == "recursive"){     
        if (missing(p0))
            return(panjer(fx = fx, x.scale = x.scale, model.freq = model.freq, echo = echo, TOL = TOL))
        else 
            return(panjer(fx, x.scale = x.scale, model.freq, p0 = p0, echo = echo, TOL = TOL))
    }
    
    if (method == "exact")
    {
        if (class(model.freq) == "numeric")
            pn <- model.freq
        else
        {
            pfreq <- match.fun(paste("p", model.freq$dist, sep = ""))
            qfreq <- match.fun(paste("q", model.freq$dist, sep = ""))
            qfreqpar <- c(p = 1 - 10e-04*TOL, model.freq$par)                
            formals(qfreq)[names(qfreqpar)] <- qfreqpar
            formals(pfreq)[names(model.freq$par)] <- model.freq$par                
            Pn <- pfreq(seq(0, qfreq()))
            pn <- c(0, diff(Pn))
        }
        return(exact(x.scale = x.scale, fx = fx, pn = pn))
    }
}
    

                     
                
            
            
       
            
            
