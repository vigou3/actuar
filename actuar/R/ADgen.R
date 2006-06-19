ADgen <- function(method = c("normal", "np2", "simulation", "recursive", "exact"),
                  model.freq, model.sev, moments, h, p0, TOL = 1e-06, ...) 
{
    if (method == "normal" | method == "np2")
    {
        fmoments <- moments
        fmoments["var"] <- c(sd = sqrt(moments["var"]))
        formals(qnorm)[names(fmoments)] <- fmoments
        
        
        x <- seq(0, qnorm(1-TOL), by = h)
        FUN <- match.fun(method)
        formals(FUN)[names(moments)] <- moments
        res <- FUN(x)
    }
    else
    {
        if (method == "simulation")
        {
           res <- simS(1/TOL, model.freq, model.sev)            
        }                                           
        else
        {
            psev <- match.fun(paste("p", model.sev$dist, sep = ""))
            qsev <- match.fun(paste("q", model.sev$dist, sep = ""))
            qsevpar <- c(p = 1 - 1e-10, model.sev$par)
            formals(qsev)[names(qsevpar)]  <- qsevpar
            formals(psev)[names(model.sev$par)] <- model.sev$par
            Fx <- psev(seq(0, qsev(), by = h))  
            fx <- c(0, diff(Fx))
            if (method == "recursive")
                ifelse(missing(p0),
                       res <- panjer(fx, model.freq$dist, as.list(model.freq$par), TOL = TOL),
                       res <- panjer(fx, model.freq$dist, as.list(model.freq$par), p0 = p0, TOL = TOL))
                
            
            if (method == "exact")
            {
                pfreq <- match.fun(paste("p", model.freq$dist, sep = ""))
                qfreq <- match.fun(paste("q", model.freq$dist, sep = ""))
                qfreqpar <- c(p = 1 - TOL, model.freq$par)                
                formals(qfreq)[names(qfreqpar)] <- qfreqpar
                formals(pfreq)[names(model.freq$par)] <- model.freq$par                
                Pn <- pfreq(seq(0, qfreq()))
                pn <- c(0, diff(Pn))                
                res <- exact(fx, pn)
            }
        }
    }
    res$h <- h
    res
}
                     
                
            
            
       
            
            
