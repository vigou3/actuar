### ===== actuar: an R package for Actuarial Science =====
###
### Bowers Gamma Approximation of the total amount of
### claims distribution
###
### See Hardy, Encyclopedia of actuarial science,
### Wiley & Sons, 2004 and Bowers, Transactions of society of
### actuaries, Expansion of probability density functions as a 
### sum of gamma densities with applications in risk theory, 1966.
###
### AUTHORS:  Vincent Goulet <vincent.goulet@act.ulaval.ca>
### and Christophe Dutang

bowergamma <- function(moments)
{
    
    ## Approximate the total amount of claims distribution using the first
    ## three or the first five (centered) moments.
    
    # moments of the standardized variable X
    alpha <- moments[1]^2 / moments[2]
    beta <- moments[1] / moments[2]
    mu3 <- beta^3 * moments[3]
#    cat("alpha", alpha, " beta ", beta,"\n")
#    cat("central mom mu1", alpha, " mu2 ", alpha, " mu3 ", mu3, "\n")
    if( length(moments) == 5)
    {
        mu4 <- beta^4 * moments[4]
        mu5 <- beta^5 * moments[5]
#        cat("mu4 ", mu4, " mu5 ", mu5, "\n")
    }
    
# coefficient of the orthogonal polynomials expansion
    A3 <- 1 / 6 * (mu3 -  2 * alpha)
    if( length(moments) == 5)
    {
        A4 <- 1 / 24 * (mu4 - 12 * mu3 - 3 * alpha^2 + 18 * alpha)
        A5 <- 1 / 120 * (mu5  - 20 * mu4 - (10 * alpha - 120) * mu3 + 60 * alpha^2 - 144 * alpha)
    }
    else
    {
        A4 <- 0
        A5 <- 0
    }
#    cat("A3 ", A3 ,"\n")
#       cat("A4 ", A4 ,"\n")
#        cat("A5 ", A5 ,"\n")
# coefficient of the Bowers gamma approximation    
    c1 <- 1 - A3 + A4 - A5
    c2 <- 3 * A3 - 4 * A4 + 5 * A5
    c3 <- -3 * A3 + 6 * A4 - 10 * A5
    c4 <- A3 - 4 * A4 + 10 * A5
    c5 <- A4 - 5 * A5
    c6 <- A5
    
    cat("c1 ", c1, " c2 ", c2, " c3 ", c3, "\n")
    cat("c4 ", c4, " c5 ",c5, " c6 ", c6, " somme ",c1+c2+c3+c4+c5+c6 ,"\n")
    
    FUN <- function(x)
    {
#print(x)
            x <- x * beta
#print(x)
            pgamma(x, alpha) * c1 + pgamma(x, alpha+1) * c2 + pgamma(x, alpha+2) * c3 
            + pgamma(x, alpha+3) * c4 + pgamma(x, alpha+4) * c5 + pgamma(x, alpha+5) * c6
    }
    
    environment(FUN) <- new.env()
    assign("mean", moments[1], envir = environment(FUN))
    assign("variance", moments[2], envir = environment(FUN))
    assign("skewness", moments[3] / sqrt(moments[2])^3, envir = environment(FUN))
    attr(FUN, "source") <- "function(x) pgamma(x, alpha) * c1 + pgamma(x, alpha+1) * c2 + pgamma(x, alpha+2) * c3 + pgamma(x, alpha+3) * c4 + pgamma(x, alpha+4) * c5 + pgamma(x, alpha+5) * c6"
    FUN
}
