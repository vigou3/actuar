aggregateDist <- function(method = c("normal", "np2", "simulation", "recursive", "exact"),
                  model.sev, model.freq, moments, x.scale = 1, n, p0, TOL = 1e-06, ...)
{
    method <- match.arg(method)
    
    if (method == "normal"){
        if (length(moments) != 2) stop("'normal' method requires the first TWO moments of the distribution")
        return(normal(moments[1], moments[2]))}
        
    if (method == "np2"){
        if (length(moments) != 3) stop("'np2' method requires the first THREE moments of the distribution")
        return(np2(moments[1], moments[2], moments[3]))}

    
    if (method == "simulation") return(simS(n, model.freq, model.sev, x.scale = x.scale))
    
    if (class(model.sev) == "numeric") fx <- model.sev
       
    else
    {
        psev <- match.fun(paste("p", model.sev$dist, sep = ""))
        qsev <- match.fun(paste("q", model.sev$dist, sep = ""))
        qsevpar <- c(p = 1 - 1e-10, model.sev$par)
        formals(qsev)[names(qsevpar)]  <- qsevpar
        formals(psev)[names(model.sev$par)] <- model.sev$par
        Fx <- psev(seq(0, qsev(), by = 1/x.scale))  
        fx <- c(0, diff(Fx))
    }
    if (method == "recursive"){     
        if (missing(p0))
            return(panjer(fx = fx, x.scale = x.scale, model.freq = model.freq, TOL = TOL))
        else 
            return(panjer(fx, x.scale = x.scale, model.freq, p0 = p0, TOL = TOL))
    }
    
    if (method == "exact")
    {
        if (class(model.freq) == "numeric")
            pn <- model.freq
        else
        {
            pfreq <- match.fun(paste("p", model.freq$dist, sep = ""))
            qfreq <- match.fun(paste("q", model.freq$dist, sep = ""))
            qfreqpar <- c(p = 1 - TOL, model.freq$par)                
            formals(qfreq)[names(qfreqpar)] <- qfreqpar
            formals(pfreq)[names(model.freq$par)] <- model.freq$par                
            Pn <- pfreq(seq(0, qfreq()))
            pn <- c(0, diff(Pn))
        }
        return(exact(x.scale = x.scale, fx = model.sev, pn = model.sev))
    }
}
    

                     
                
            
            
       
            
            
