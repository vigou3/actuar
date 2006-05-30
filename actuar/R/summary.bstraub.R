summary.bstraub <- function(object, ...)
    structure(object, class = c("summary.bstraub", class(object)))

print.summary.bstraub <- function(x, ...)
{
    cat("\nCredibility model:", x$model, "\n\n")
    cat("Structure Parameters Estimators\n\n")
    cat("  Collective premium:        ", x$collective, "\n")
    cat("  Average contract variance: ", x$s2,"\n")
    cat("  Portfolio heterogeneity:   ", x$unbiased, " (unbiased)\n")
    cat("                             ",
        ifelse(is.null(x$iterative), "NULL     ", x$iterative),
        " (iterative)\n")
    cat("  Credibility constant:      ",
        x$s2 / ifelse(is.null(x$iterative), x$unbiased, x$iterative),
        "\n\n")
    cat("Detailed premiums\n\n")
    cred <- cbind(1:x$ncontracts, x$individual, x$weights, x$cred, x$premium)
    colnames(cred) <- c(" Contract", "Ind. premium", "Weight",
                        "Cred. factor", "Cred. premium")
    rownames(cred) <- rep("", x$ncontracts)
    print(cred)
    invisible(x)
}
