### ===== actuar: an R package for Actuarial Science =====
### 
### A summary method for class 'simpf'
###
### AUTHORS:  Louis-Philippe Pouliot, Vincent Goulet <vincent.goulet@act.ulaval.ca>

summary.simpf <- function(object,
                          contracts = seq(length = dim(object$data)[1]),
                          years = seq(length = dim(object$data)[2]),
                          ...)
{
    object$contracts <- contracts
    object$years <- years
    class(object) <- "summary.simpf"
    object    
}

print.summary.simpf <- function(x, contracts = x$contracts, years = x$years, ...)
{
    cat("\nIndividual claim amounts\n")
    for (i in contracts)
    {
        cat("\n Contract",i,"\n")
        sapply(paste("   Year",
                 years,
                 ":",
                 format(x$data[i,years])),
           cat,
           fill = TRUE)
    }
    invisible(x)
}

