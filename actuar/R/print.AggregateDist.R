print.AggregateDist <- function(x, ...)
{
    cat("\nEmpirical CDF for the aggregate claims\n\n")
    cat("Call:\n")
    cat("  ",deparse(x$call), "\n\n", fill = TRUE)
    cat("Data:  (", length(x$Fs), "obs. )\n\n")

    digits <- 4
    numform <- function(x) paste(formatC(x, dig = digits, width = 5), collapse = ", ")
    n <- length(x$Fs)
    i1 <- 1:min(3, n)
    i2 <- if (n >= 4) 
        max(4, n - 1):n
    else integer(0)
    xx <- 1:n


    
    cat("    S:  ", numform(xx[i1]), if (n > 3) 
        ", ", if (n > 5) 
        " ..., ", numform(xx[i2]), "\n", sep = "")

    cat("  DDF: ", numform(x$fs[i1]), if (n > 3) 
        ", ", if (n > 5) 
        " ..., ", numform(x$fs[i2]), "\n", sep = "")

    cat("  CDF: ", numform(x$Fs[i1]), if (n > 3) 
        ", ", if (n > 5) 
        " ..., ", numform(x$Fs[i2]), "\n", sep = "")
    cat("\n\n")

    
    cat("Summary\n")
    q <- quantile.AggregateDist(x, p = c(0.25, 0.5, 0.75))
    expectation <- mean.AggregateDist(x)
    summary <- c(0, q[1:2], expectation, q[3], n)
    names(summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    print(summary)
    res <- list(FUN = approxfun(x$Fs))
    res
}


    
