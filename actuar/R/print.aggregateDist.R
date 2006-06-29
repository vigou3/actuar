print.aggregateDist <- function(x, ...)
{
    cat("\nAggregate Claim Amounts Empirical CDF\n")
    cat("  ", label <- get("label", environment(x)), "\n\n")
    #cat("  ", label <- comment(x), "\n\n")
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


    #cat("      S:  ", numform(xx[i1]), if (n > 3) 
    #    ", ", if (n > 5) 
    #   " ..., ", numform(xx[i2]), "\n", sep = "")

    #cat("PMF/PDF: ", numform(x$fs[i1]), if (n > 3) 
    #    ", ", if (n > 5) 
    #    " ..., ", numform(x$fs[i2]), "\n", sep = "")

    #cat("    CDF: ", numform(x$Fs[i1]), if (n > 3) 
    #    ", ", if (n > 5) 
    #    " ..., ", numform(x$Fs[i2]), "\n", sep = "")
    #cat("\n")

    
    #cat("Summary\n")
    #q <- quantile.aggregateDist(x, p = c(0.25, 0.5, 0.75))
    #expectation <- mean.aggregateDist(x)
    #summary <- c(0, q[1:2], expectation, q[3], n)
    #names(summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    #print(summary)

    
