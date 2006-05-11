"print.bstraub" <-
function(x,...)
{
    if (x$unspec.weights)
    {
        model <- "Bühlmann"
        underline <- "-----------------------------------------"
    }
    else
    {
        model <- "Bühlmann-Straub"
        underline <- "------------------------------------------------"
    }
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Credibility Premiums using", model, "model \n")
    cat(underline,"\n")

    sapply(paste("Contract ",
                 format(1:x$ncontracts, justify="right"),
                 ": ",
                 format(x$premiums),
                 sep = ""),
           cat, fill = TRUE)   
    invisible(x)
}

