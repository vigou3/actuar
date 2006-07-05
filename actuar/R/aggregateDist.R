### ===== actuar: an R package for Actuarial Science =====
###
### Use one of five methods to compute the aggregate claims
### distribution of a portfolio over one year given a frequency and a
### severity model or the true moments of the distribution.
###
### AUTHORS:  Louis-Philippe Pouliot, Vincent Goulet <vincent.goulet@act.ulaval.ca>


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

print.aggregateDist <- function(x, ...)
{
    cat("\nAggregate Claim Amounts Empirical CDF\n")
    cat("  ", label <- comment(x), "\n\n")
    if (label %in% c("Direct calculation", "Recursive method approximation"))
        cat("Discretization step :", get("x.scale", envir = environment(x)), "\n\n")
    
    cat("Call:\n")
    print(get("call", envir = environment(x)))
    cat("\n")

    if (label %in% c("Direct calculation",
                      "Recursive method approximation",
                      "Approximation by simulation"))
    {
        n <- length(get("x", environment(x)))
        cat("Data:  (", n, "obs. )\n")
        numform <- function(x) paste(formatC(x, dig = 4, width = 5), collapse = ", ")
        i1 <- 1:min(3, n)
        i2 <- if (n >= 4)
            max(4, n - 1):n
        else integer(0)
        xx <- eval(expression(x), env = environment(x))
        cat(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3) 
        ", ", if (n > 5) 
        " ..., ", numform(xx[i2]), "\n", sep = "")
        cat("\n")
    } 
    if (label %in% c("Normal approximation",
                      "Normal Power approximation"))
    {
        cat(attr(x, "source"), "\n")
    }
    print(environment(x))
    cat("Class attribute:\n")
    print(attr(x, "class"))   
}

plot.aggregateDist <- function(x, xlim, ...)
{
    ## Function plot() is used for the discretized
    ## CDFs and function curve() in the continuous cases.
    main <- "Aggregate Claims Distribution"
    
    if ("stepfun" %in% class(x)){

        ## Method for class 'ecdf' will most
        ## probably be used.
        NextMethod(ylab = "", main = main) 
    }
    else
    {
        ## Limits for the x-axis are supplied if none are given
        ## in argument.
        
        if (missing(xlim)){
            mean <- get("mean", environment(x))
            sd <- sqrt(get("var", environment(x)))
            xlim <- c(mean - 3*sd, mean + 3*sd)
        }
        curve(x, main = main, ylab = "", xlim = xlim, ylim = c(0,1))
    }
    mtext(expression(F*scriptstyle(s)(x)), side = 2, line = 2)
    mtext(comment(x), line = 0)    
}

summary.aggregateDist <- function(object, ...)
    {structure(object, class = c("summary.aggregateDist", class(object)))}

print.summary.aggregateDist <- function(x, ...)
{
    cat(ifelse(comment(x) %in% c("Normal approximation","Normal Power approximation"),
               "Aggregate Claim Amounts CDF:\n",
               "Aggregate Claim Amounts Empirical CDF:\n"))
    q <- quantile(x, p = c(0.25, 0.5, 0.75))
    expectation <- mean(x)
    
    if (comment(x) %in% c("Normal approximation","Normal Power approximation"))
    {
        min <- 0; max <- NA
    }
    else
    {
        max <- tail(eval(expression(x), environment(x)), 1)
        min <- head(eval(expression(x), environment(x)), 1)
    }
    res <- c(min, q[1:2], expectation, q[3], max)
    names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    print(res)
}

mean.aggregateDist <- function(x, ...)
{
    label <- comment(x)	

    ## Simply return the value of the true mean
    ## given in argument in the case of the Normal
    ## and Normal Power approximations.
    
    if (label %in% c("Normal approximation",
                      "Normal Power approximation"))
        return(eval(expression(mean), environment(x)))
    
    else
        return(crossprod(get("x", environment(x)),
                         c(0, diff(eval(expression(y), environment(x))))))
}
        

    

                     
                
            
            
       
            
            
