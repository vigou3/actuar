simAD <- function(method, model.freq, model.sev, moments, h, p0, TOL = 1e-06, ...) ###que faire avec p0
{
    if (method == "normal" | method == "np2")
    {
        lim = qnorm(1-TOL)
        x <- seq(-lim, lim, by = h)
        par <- c(x = 1/h, moments)
        FUN <- match.fun(method)
        formals(FUN) <- par
        res <- FUN
    }
    else
    {
        if (method == "sim")
        {
           X <- simS(1/TOL, model.freq, model.sev) ### décider quoi faire avec la distribution
       }                                          ### appeler density()???
        else
        {
            dsev <- match.fun(paste("d", model.sev$dist, sep = ""))
            qsev <- match.fun(paste("q", model.sev$dist, sep = ""))
            qsevpar <- c(p = 1 - TOL, model.sev$par)
            formals(qsev)[names(qsevpar)]  <- qsevpar
            formals(dsev)[names(model.sev$par)] <- model.sev$par
            fx <- dsev(seq(0, qsev(), by = h))
            
            if (method == "recursive")
                res <- panjer(fx, model.freq$dist, model.freq$par, p0, TOL)
            
            if (method == "exact")
            {
                dfreq <- match.fun(paste("d", model.freq$dist, sep = ""))
                qfreq <- match.fun(paste("q", model.freq$dist, sep = ""))
                qfreqpar <- c(p = 1 - TOL, model.freq$par)
                formals(qfreq)[names(qfreqpar)] <- qfreqpar
                formals(dfreq)[names(model.freq$par)] <- model.freq$par
                pn <- dfreq(seq(0, qfreq()))
                res <- exact(fx, pn)
            }
        }
    }
}
            
                
                
            
            
       
            
            
